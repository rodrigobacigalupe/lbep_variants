#! /usr/bin/python
#
# variant_calling.py
#
# Author: Rodrigo Bacigalupe (built from Paul McAdam's script)
# This file is part of the LBEP bacterial genome mapping software.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

##################
# IMPORT MODULES #
##################

import os, sys, datetime
from optparse import OptionParser, OptionGroup
from random import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#########################
# SET PROGRAM LOCATIONS #
#########################

BWA_DIR="path_to_software/bwa/bwa-0.5.9/"
SAMTOOLS_DIR="path_to_software/samtools/samtools-0.1.18/"
BAMTOOLS_DIR="path_to_software/bamtools/"
BCFTOOLS_DIR="path_to_software/samtools/samtools-0.1.18/bcftools/"
PICARDTOOLS_DIR="path_to_software/picard/picard-tools-1.90/"
ADDITIONAL_SCRIPTS_DIR="path_to_scripts/"
JAVA_DIR="path_to_software/jre1.7.0_40/bin/"
SNPEFF_DIR="path_to_software/snpEff/"
GTAK_DIR="path_to_software/GenomeAnalysisTK-2.7-4/"
PINDEL_DIR="path_to_software/pindel/"

#################
# ERROR MESSAGE #
#################

def errorMessage(error):
        print "\n**ERROR: ", error
        print "\nFor help use -h or --help\n"
        sys.exit()

###########################
# USER PROVIDED ARGUMENTS #
###########################

# Parse options for running programs

def user_options():
        usage = "%prog [options] <directory of fastq files>"
        parser = OptionParser(usage=usage)
        #do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
        parser.disable_interspersed_args()

        group = OptionGroup(parser, "Required options")
        group.add_option("-n", "--project", action="store", dest="projectname", help="Name for project", default="")
        group.add_option("-r", "--reference", action="store", dest="reference", help="Reference genome in fasta format", default="")
        group.add_option("-b", "--genbank", action="store", dest="genbank", help="Genbank file", default="")
        parser.add_option_group(group)

        group = OptionGroup(parser, "Mapping options")
        group.add_option("-q", "--quality", action="store", type="int", dest="quality", help="Minimum mapping quality [default= %default]", default=30, metavar="INT")
        group.add_option("-s", "--swa", action="store", type="choice", dest="smithwaterman", choices=["True", "False"], help="enable/disable Smith-Waterman algorithm [default= %default]", default="True")
        group.add_option("-d", "--depth", action="store", type="int", dest="depth", help="Minimum depth to call base [default= %default]", default=3, metavar="INT")
        group.add_option("-p", "--percent", action="store", type="float", dest="percent", help="Minimum percent of reads containing snp required to call snp [default= %default]", default=0.66)
        parser.add_option_group(group)

        group = OptionGroup(parser, "GATK options")
        group.add_option("-c", "--calls", action="store", dest="calls", help="Minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called. [Default= %default]", default="30", metavar="INT")
        group.add_option("-e", "--emit", action="store", dest="emit", help="Minimum confidence threshold (phred-scaled) at which the program should emit sites that appear to be possibly variant. [Default= %default]", default=10, metavar=
"INT")
        parser.add_option_group(group)

        group = OptionGroup(parser, "Large deletions options")
        group.add_option("-g", "--genome", action="store", dest="genome", help="Genome file for BEDtools containing two columns, one with the name of the genome and another one with its size.", default="")
        #group.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
        parser.add_option_group(group)

        group = OptionGroup(parser, "Misc. options")
        group.add_option("-t", "--threads", type='int', action="store", dest="threads", help="number of threads to use", default=1, metavar="INT")
        parser.add_option_group(group)

        return parser.parse_args()

########################
# CHECK USER ARGUMENTS #
########################

def check_input_validity(options, args):
        if options.reference=="":
                errorMessage('No reference file provided (-r)')
        if options.projectname=="":
                errorMessage('No project name provided (-n)')
        if options.genome=="":
                errorMessage('No genome file provided (-g)')
        elif options.quality>60:
                errorMessage('Mapping quality (-q) set unrealistically high (>60)')
        elif options.threads > 5:
                errorMessage('Too many threads selected, please consider other users (>5)')

        return

#####################
# MAPPING FUNCTIONS #
#####################

def returnStrains(cwd):
        STRAIN_DIRECTORY=cwd
        strain_folders=os.listdir(STRAIN_DIRECTORY)
        strain_folders.sort()
        strain_folders1=[]
        for strain in strain_folders:
                if strain.startswith('strain'):
                        strain_folders1.append(strain)
        return strain_folders

class fastqMapping:
        def __init__(self, fastq='', name='', number=0):
                self.fastq=fastq
                self.name=name
                self.runname=''
                self.fastqdir=''
                self.number=number

        def runBWA(self, reference):
                print "\nRunning BWA aln on "+self.name         #Map reads to genome
                print "\nselffastqdir "+self.fastqdir   #Map reads to genome
                print "\nself.runname"
                os.system(BWA_DIR+"bwa aln "+options.reference+' '+self.fastqdir+self.name+'R1.fastq.gz > '+self.runname+'/tmp_1.sai')
                os.system(BWA_DIR+'bwa aln '+options.reference+' '+self.fastqdir+self.name+'R2.fastq.gz > '+self.runname+'/tmp_2.sai')

                print "\nRunning BWA sampe on "+self.name       #Combine and produce SAM file
                if options.smithwaterman=="True":
                        os.system(BWA_DIR+'bwa sampe '+reference+' '+self.runname+'/tmp_1.sai '+self.runname+'/tmp_2.sai '+self.fastqdir+self.name+'R1.fastq.gz '+self.fastqdir+self.name+'R2.fastq.gz > '+self.runname+'/tmp.sam')
                if options.smithwaterman=="False":
                        os.system(BWA_DIR+'bwa sampe -s '+reference+' '+self.runname+'/tmp_1.sai '+self.runname+'/tmp_2.sai '+self.fastqdir+self.name+'R1.fastq.gz '+self.fastqdir+self.name+'R2.fastq.gz > '+self.runname+'/tmp.sam')

                print "\nProducing BAM file for "+self.name     #produce the BAM file and filter quality
                os.system(JAVA_DIR+'java -jar '+PICARDTOOLS_DIR+'SamFormatConverter.jar I='+self.runname+'/tmp.sam '+'O='+self.runname+'/tmp.bam VALIDATION_STRINGENCY=SILENT')
                os.system(BAMTOOLS_DIR+'bin/bamtools filter -in '+self.runname+'/tmp.bam -out '+self.runname+'/tmpout.bam -script '+BAMTOOLS_DIR+'scripts/filter.json')
                #os.system(SAMTOOLS_DIR+'samtools view -bt1 -q '+str(options.quality)+' '+self.runname+"/tmp.sam > "+self.runname+"/tmp.bam")

                return

        def bamFormatting(self, reference):
                print 'Adding read group for '+self.name        #index the BAM file
                os.system(JAVA_DIR+'java -jar '+PICARDTOOLS_DIR+'AddOrReplaceReadGroups.jar INPUT='+self.runname+'/tmpout.bam OUTPUT='+self.runname+'/tmprg.bam LB=1 PL=1 PU=1 SM=1')
                print 'Sorting BAM file for '+self.name #sort by coordinate Using Piccard
                os.system(JAVA_DIR+'java -Xms1024m -Xmx1024m -jar '+PICARDTOOLS_DIR+'SortSam.jar INPUT='+self.runname+'/tmprg.bam OUTPUT='+self.runname+'/tmprgsort.bam SORT_ORDER=coordinate')
                print 'Marking duplicates '+self.name           #Mark duplicates using Piccard
                os.system(JAVA_DIR+'java -jar '+PICARDTOOLS_DIR+'MarkDuplicates.jar INPUT='+self.runname+'/tmprgsort.bam OUTPUT='+self.runname+'/dedup_reads.bam METRICS_FILE='+self.runname+'/metrics.txt VALIDATION_STRINGENCY=SILENT ')
                os.system('rm '+self.runname+'/tmp*')
                print 'Indexing BAM file for '+self.name                #sort the BAM file
                os.system(JAVA_DIR+'java -jar '+PICARDTOOLS_DIR+'BuildBamIndex.jar INPUT='+self.runname+'/dedup_reads.bam')

                return

        def runGATK(self, reference):
                #Perform local realignment around indels
                #First create a target list of intervals to be realigned
                print 'Performing local realignment around indels for '+self.name
                print 'Creating a target list of intervals to be realigned'
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T RealignerTargetCreator -R '+options.reference+' -I '+self.runname+'/dedup_reads.bam -o '+self.runname+'/target_intervals.list')
                #Secondly perform realignment of the target intervals
                print 'Performing realignment of the target intervals'
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T IndelRealigner -R '+options.reference+' -I '+self.runname+'/dedup_reads.bam -targetIntervals '+self.runname+'/target_intervals.list'+' -o '+self.runname+'/
realigned_reads.bam')
                #Index the new bam file and remove dedup
                print 'Indexing BAM file for '+self.name
                os.system(JAVA_DIR+'java -jar '+PICARDTOOLS_DIR+'BuildBamIndex.jar INPUT='+self.runname+'/realigned_reads.bam')
                os.system('rm '+self.runname+'/dedup_reads.*')

                #Variant discovery with UnifiedGenotyper
                print 'Calling variants for '+self.name
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T UnifiedGenotyper -R '+options.reference+' -I '+self.runname+'/realigned_reads.bam -rf BadCigar -ploidy 1 -glm BOTH -stand_call_conf '+str(options.calls)+'
-stand_emit_conf '+str(options.emit)+' -o '+self.runname+'/raw_haploid_variants.vcf')

                #Variant recalibration (applying hard filters to a call set)
                print 'Applying hard filters to '+self.name
                print 'Extracting SNPs from the call set of '+self.name
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T SelectVariants -R '+options.reference+' -V '+self.runname+'/raw_haploid_variants.vcf -selectType SNP -o '+self.runname+'/raw_snps.vcf')
                print 'Applying filters to the SNP call set for '+self.name
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T VariantFiltration -R '+options.reference+' -V '+self.runname+'/raw_snps.vcf --clusterSize  3 --clusterWindowSize  10 --filterExpression "QD < 2.0 || FS > 6
0.0 || MQ < 40.0 || QUAL <30" -filterName "my_snp_filter"  -o '+self.runname+'/'+self.runname+'filtered_snps.vcf')
                print 'Extracting indels from the call set of '+self.name
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T SelectVariants -R '+options.reference+' -V '+self.runname+'/raw_haploid_variants.vcf -selectType INDEL -o '+self.runname+'/raw_indels.vcf')
                print 'Applying filters to the indels call set for '+self.name
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T VariantFiltration -R '+options.reference+' -V '+self.runname+'/raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0" -filterName "my_indel_filter"  -o
 '+self.runname+'/'+self.runname+'filtered_indels.vcf')
                print 'Merge SNPs and indels files into variant files'
                os.system(JAVA_DIR+'java -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T CombineVariants -R '+options.reference+' --variant '+self.runname+'/'+self.runname+'filtered_snps.vcf --variant '+self.runname+'/'+self.runname+'filtered_i
ndels.vcf -o '+self.runname+'/'+self.runname+'filtered_variants.vcf -genotypeMergeOptions UNIQUIFY')
                #Variant discovery
                print 'Producing consensus fasta for '+self.name
                os.system(JAVA_DIR+'java -Xmx2g -jar '+GTAK_DIR+'GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R '+options.reference+' -o '+self.runname+"/"+self.name+'consensus.fa --variant '+self.runname+'/'+self.runname+'filte
red_variants.vcf')
                os.system('mv '+self.runname+"/"+self.name+'consensus.fa '+consensus_directory)
                os.system('mv '+self.runname+"/"+self.name+'filtered_*.vcf '+gatk_directory)
                return

        def runsamtools(self, reference):
                #Produces vcf file and consensus sequence
                print 'Producing vcf for '+self.name
                os.system(SAMTOOLS_DIR+'samtools mpileup -C 50 -uf '+options.reference+' '+self.runname+'/realigned_reads.bam | '+BCFTOOLS_DIR+'bcftools view -cg - > '+self.runname+"/"+self.name+'.vcf')
                print 'Producing consensus fasta for '+self.name
                os.system(ADDITIONAL_SCRIPTS_DIR+'vcf_to_consensus_sequence.py -r '+options.reference+' -d '+str(options.depth)+' -p '+str(options.percent)+' '+self.runname+"/"+self.name+'.vcf')
                os.system('mv '+self.runname+"/"+self.name+'.consensus.fa '+consensus_directory+self.name+'consensus_samtools.fa')
                #os.system('mv '+self.runname+"/"+self.name+'.vcf '+vcf_directory+'consensus_vcf') #This step is temporarily commented.
                os.system('rm '+self.runname+"/"+self.name+'.vcf')                                 #You can comment this one instead of the previous one.
                os.system('mv '+self.runname+"/"+self.name+'_snp.vcf '+samtools_directory)
                return

        def runpindel(self, reference):
                print 'Running pindel for '+self.name   #Run pindel
                print 'Extracting useful reads for Pindel from BAM files'
                os.system('perl '+PINDEL_DIR+'bam2pindel.pl -i '+intermediate_files+self.runname+'/realigned_reads.bam'+' -o '+self.runname+' -s '+self.runname+' -om -pi 250')
                print 'Calling indels and structural variations with Pindel '+self.name
                os.system(PINDEL_DIR+'pindel -f '+options.reference+' -p '+self.runname+'* -c ALL -o '+self.runname)
                print 'Filtering files... '             #Removing empty files
                os.system('find ./ -name "strain*" -size 0 -exec rm {} \;')
                reference_name = str(options.reference.split('/')[-1]).split('.')[0]
                date_time = str(datetime.date.today()).replace("-", "")
                os.system('for i in '+self.runname+'_*; do    '+PINDEL_DIR+'pindel2vcf -p $i -r '+options.reference+' -R '+reference_name+' -d '+date_time+' -G; done')
                os.system('grep "ChrID" '+self.runname+'_* | cut -f 1,2,3,5,6,9 > '+self.runname)
                os.system('rm *.txt.vcf *__SI *__D')

                return

        def largedeletions(self, reference):
                print 'Identifying large deletions for '+self.name      #Identification of large deletions
                print 'Converting bam files to bed files for '+options.reference
                os.system('bedtools bamtobed -i '+intermediate_files+self.runname+'/realigned_reads.bam > '+self.runname+'tmp_reads.bed')
                print 'Calculate the number of reads and coverage density for '+self.name
                os.system('coverageBed -a '+self.runname+'tmp_reads.bed'+' -b ../1000bpwindow.bed | sortBed -i stdin > '+self.runname+'.tmp.windows.bed.coverage;')
                print 'Filter the file to identify low coverages'
                os.system("awk '$7<=0.5' "+self.runname+'.tmp.windows.bed.coverage > '+self.runname+'low_coverage')
                os.system('rm *tmp*')
                os.system(ADDITIONAL_SCRIPTS_DIR+'large_deletions.py '+self.runname+'low_coverage '+options.genbank)
                return

        def largeinsertions(self, reference):
                print 'Identifying large insertions for '+self.name     #Identification of large insertions
                print 'Searching for unmapped reads for '+options.reference
                os.system(SAMTOOLS_DIR+'samtools view -F 2 '+intermediate_files+self.runname+'/realigned_reads.bam | awk \'{if ($3 == "*" ) print ">"$1"\n"$10}\' > '+self.runname+'unmapped_reads.fa')
                print 'Calculate the number of reads and coverage density for '+self.name
                #os.system('coverageBed -a '+self.runname+'tmp_reads.bed'+' -b ../1000bpwindow.bed | sortBed -i stdin > '+self.runname+'.tmp.windows.bed.coverage;')
                #print 'Filter the file to identify low coverages'
                #os.system("awk '$7<=0.5' "+self.runname+'.tmp.windows.bed.coverage > '+self.runname+'low_coverage')
                #os.system('grep "ChrID" '+self.runname+'_* | cut -f 1,2,3,5,6,9 > '+self.runname)
                #os.system('rm *.txt *__*')
                return

def runSNPeff(reference, vcf):
        strain_snpeff_dir=snpeff_directory+strain.name+'/'
        os.mkdir(strain_snpeff_dir)
        os.system(SNPEFF_DIR+'snpEff.jar eff -c '+SNPEFF_DIR+'snpEff.config -o txt -ud 0 -a 50 -s '+strain_snpeff_dir+strain.name+'.html '+reference+' '+vcf+' > '+strain_snpeff_dir+strain.name+'.snpeff')
        return

def makeSnpSpreadsheet(reference, projectname):
        os.mkdir('./temp')
        os.chdir('./temp')
        snpeffs=os.listdir(snpeff_directory)
        for snpeff in snpeffs:
                os.system('cp '+snpeff_directory+snpeff+'/*.snpeff ./')
        os.system(ADDITIONAL_SCRIPTS_DIR+'snp_spreadsheet_maker.py -r '+reference+' -n '+projectname+' ./')
        os.system('mv '+projectname+'.xls ../'+projectname+'_snps.xls')
        os.chdir('../')
        os.system('rm -r ./temp')
        return


########
# MAIN #
########

if __name__ == '__main__':
        # Parse options and check input
        (options, args)=user_options()
        check_input_validity(options, args)

        # Set up directories and create folders
        os.chdir(args[0])
        cwd=os.getcwd()
        intermediate_files=cwd+'/Intermediate_files_'+options.projectname+'/'
        consensus_directory=cwd+'/Consensus_'+options.projectname+'/'
        samtools_directory=cwd+'/Samtools_'+options.projectname+'/'
        gatk_directory=cwd+'/GATK_'+options.projectname+'/'
        vcf_directory=cwd+'/VCFs_'+options.projectname+'/'
        snpeff_directory=cwd+'/snpEff_'+options.projectname+'/'
        pindel=cwd+'/Pindel_'+options.projectname+'/'
        large_deletions=cwd+'/Large_deletions_'+options.projectname+'/'
        large_insertions=cwd+'/Large_insertions_'+options.projectname+'/'

        folders = returnStrains(cwd)

        os.mkdir(intermediate_files)
        os.mkdir(consensus_directory)
        os.mkdir(samtools_directory)
        os.mkdir(gatk_directory)
        os.mkdir(vcf_directory)
        os.mkdir(vcf_directory+'core_snp_vcf')
        #os.mkdir(vcf_directory+'snp_vcf')
        os.mkdir(vcf_directory+'consensus_vcf')
        os.mkdir(snpeff_directory)
        os.mkdir(pindel)
        os.mkdir(large_deletions)
        os.mkdir(large_insertions)

        #print 'Loading required modules'
        #os.system('module add /opt/modules/modulefiles/apps/gcc/BEDTools/2.17.0')
        print '\nValidating input files...'

        strains=[]

        count=0

        for folder in folders:
                if os.path.isdir(folder):
                        files=os.listdir(folder)
                        file_count=len(files)
                        if file_count!=2:               #Check there are two files in each folder, one forward and one reverse read file
                                print '**WARNING: Incorrect number of files in '+folder
                                sys.exit()
                        for f in files:
#                               if f[-3:]!='.fq':       #Check files are fastq, must have .fq suffix
#                                       print '**WARNING: Files are not fastq in '+folder
#                                       sys.exit()
                                assembly_name=f[:-11]

                        if not os.path.isdir(intermediate_files+assembly_name):
                                print 'Creating '+intermediate_files+assembly_name
                                os.mkdir(intermediate_files+assembly_name)
                                os.mkdir(large_insertions+assembly_name)

                        strains.append(fastqMapping())
                        strains[count].number=str(count+1)
                        strains[count].runname=assembly_name
                        strains[count].name=assembly_name
                        strains[count].fastqdir=cwd+'/'+folder+'/'

                        count+=1
        print str(len(strains))

        if len(strains)==0:
                print '\n**ERROR: No valid input files'
                sys.exit()

        """
        # Prepare a reference for use with BWA and GATK
        # Generate the BWA index
        os.system(BWA_DIR+'bwa index -a is'+options.reference)
        # Generate the fasta file index
        os.system(SAMTOOLS_DIR+'samtools faidx '+options.reference)
        # Generate the sequence dictionary
        os.system(JAVA_DIR+'java -jar '+PICARDTOOLS_DIR+'CreateSequenceDictionary.jar REFERENCE='+options.reference+' OUTPUT='+options.reference+'.dict')
        """

        #Create a window of 1000 bp for identifying large deletions
        os.system('bedtools makewindows -g '+options.genome+' -w 1000 > 1000bpwindow.bed')

        for strain in strains:
                os.chdir(intermediate_files)
                strain.runBWA(options.reference)
                strain.bamFormatting(options.reference)
                strain.runGATK(options.reference)
                strain.runsamtools(options.reference)
                os.chdir(pindel)
                strain.runpindel(options.reference)
                os.chdir(large_deletions)
                strain.largedeletions(options.reference)
                #os.chdir(large_insertions)
                #strain.largeinsertions(options.reference)
                os.chdir(cwd)

                #coreVcfSnps(core_position_set, vcf_directory+'snp_vcf/'+strain.name)
                #runSNPeff(options.reference.split('/')[-1].split('.')[0], vcf_directory+'core_snp_vcf/'+strain.name+'_core_snp.vcf')

        #makeSnpSpreadsheet(options.reference, options.projectname)
        os.system('rm 1000bpwindow.bed')
        os.mkdir(options.projectname+'_results/')
        os.system('mv *'+options.projectname+' '+options.projectname+'_results/')

