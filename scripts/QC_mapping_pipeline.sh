#Usage: bash QC_mapping_pipeline.sh -n basename_of_each_library_fastq_file -g path/to/genome.fasta
#The pipeline doesn't need to be modified if your basenames fit the format of XXXX_R1.fastq or XXXX_R2.fastq (where the XXXX can be anything as long as it doesn't contain spaces or underscores, the only underscore should be before the R1 or R2).
#If you have a different file format, then you need to modify the MSI job submission script (placed in the analysis directory).
#Use the rename.py script in the scripts directory to rename files from UMGC to fit this format.
#Make sure your referenced genome in the reference_genomes folder is indexed with BWA and samtools. Here are the commands: bwa index reference.fasta and samtools faidx reference.fasta
#Read through script carefully before running to ensure all assumptions and settings match your system.

#Make a project directory (named whatever you like), and then set up the following subdirectory structures for this pipeline:
#./raw_data
#./reference_genomes
#./scripts     Make sure to add this directory to your $PATH
#./analysis     Place the MSI job script and submit it from inside this directory
 
#Check that the correct parameters were given, and if so assign the variables of LIBRARY and GENOME
while [[ $# -gt 1 ]]
do
	key="$1"
	
	case $key in
		-n|--NameofFile) #Name of the input sequencing library file prefix
		LIBRARY="$2"
		shift
		;;
		-g|--genome) #Give the path to the indexed genome fasta file
		GENOME="$2"
		shift
		;;
		#Exit here if -n, -g, or -l options were not given correctly
		*) echo "Usage: SNP_pipeline.sh -n path/to/filename_prefix -g path/to/genome.fasta. Unknown option: $1 " >&2 ; exit 1
		;;
	esac
	shift
done

#Give error and exit if the variables were not set correctly
if [ -z $LIBRARY ]; then echo "No library directory given! Usage: QC_mapping_pipeline.sh -n LIBname -g Path/to/Genome"; exit 1; fi
if [ -z $GENOME ]; then echo "No reference genome given! Usage: QC_mapping_pipeline.sh -n LIBname -g Path/to/Genome"; exit 1 ; fi

#Print the values assigned to the variables
echo "LIBRARY" = "${LIBRARY}"
echo "GENOME" = "${GENOME}"

#Identify the raw read files. 
fasR=$(ls $LIBRARY_R1.fastq)
fasL=$(ls $LIBRARY_R2.fastq)

#If input files are not found then exit.
if [ ! -f $fasR ] || [ ! -f $fasL ]
then
	echo "Reads from "${LIBRARY}" could not be found"; exit 1
fi


####################################
######       SNP calling      ###### 
####################################


#Perform filtering and trimming of the raw reads, followed by QC
                
mkdir ./trimmed_reads #output folder for the filtered and trimmed reads

#$TRIMMOMATIC is an MSI specific environment variable. On your system, just replace with the path to the Trimmomatic executables.
#Change the threads to match your system
#Change the adapter and other parameters to fit your requirements

java -jar $TRIMMOMATIC/trimmomatic.jar PE -threads 24 -trimlog ./trimmed_reads/trim_log $fasR $fasL ./trimmed_reads/${LIBRARY}_trim_R1.PE ./trimmed_reads/${LIBRARY}_trim_R1.SE ./trimmed_reads/${LIBRARY}_trim_R2.PE ./trimmed_reads/${LIBRARY}_trim_R2.SE ILLUMINACLIP:$TRIMMOMATIC/adapters/all_illumina_adapters.fa:2:30:10

#Make output dirs for results of fastx toolkit
mkdir -p ./trimmed_reads/{statistics,plots/quality_plots,plots/nt_distribution}         

#There should be very few SE reads since our trimming was not very aggressive.
#Whatever SE files are not used into mapping.
fastx_quality_stats -i ./trimmed_reads/$LIBRARY_trim_R1.PE -o ./trimmed_reads/statistics/$LIBRARY_R1_PE_stats.txt
fastx_quality_stats -i ./trimmed_reads/$LIBRARY_trim_R1.SE -o ./trimmed_reads/statistics/$LIBRARY_R1_SE_stats.txt
fastx_quality_stats -i ./trimmed_reads/$LIBRARY_trim_R2.PE -o ./trimmed_reads/statistics/$LIBRARY_R2_PE_stats.txt
fastx_quality_stats -i ./trimmed_reads/$LIBRARY_trim_R2.SE -o ./trimmed_reads/statistics/$LIBRARY_R2_SE_stats.txt

fastq_quality_boxplot_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R1_PE_stats.txt -o ./trimmed_reads/plots/quality_plots/$LIBRARY_R1_PE_quality.png
fastq_quality_boxplot_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R1_SE_stats.txt -o ./trimmed_reads/plots/quality_plots/$LIBRARY_R1_SE_quality.png
fastq_quality_boxplot_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R2_PE_stats.txt -o ./trimmed_reads/plots/quality_plots/$LIBRARY_R2_PE_quality.png
fastq_quality_boxplot_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R2_SE_stats.txt -o ./trimmed_reads/plots/quality_plots/$LIBRARY_R2_SE_quality.png
                
fastx_nucleotide_distribution_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R1_PE_stats.txt -o ./trimmed_reads/plots/nt_distribution/$LIBRARY_R1_PE_nt_distr.png
fastx_nucleotide_distribution_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R1_SE_stats.txt -o ./trimmed_reads/plots/nt_distribution/$LIBRARY_R1_SE_nt_distr.png
fastx_nucleotide_distribution_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R2_PE_stats.txt -o ./trimmed_reads/plots/nt_distribution/$LIBRARY_R2_PE_nt_distr.png
fastx_nucleotide_distribution_graph.sh -i ./trimmed_reads/statistics/$LIBRARY_R2_SE_stats.txt -o ./trimmed_reads/plots/nt_distribution/$LIBRARY_R2_SE_nt_distr.png

#Align against the reference genome (previously indexed), and then convert to bam, remove duplicates, sort and index.
mkdir ./mapped_reads #output directory for mapped reads

#Assign unique readgroups with the -R flag so that samples can be used in pooled SNP calling later
bwa mem -t 24 -R "@RG\tID:${LIBRARY}\tSM:${LIBRARY}" $GENOME ./trimmed_reads/$LIBRARY_trim_R1.PE ./trimmed_reads/$LIBRARY_trim_R2.PE | samtools view -bh -@ 24 - | samtools rmdup - | samtools sort -o ./mapped_reads/$LIBRARY.rm.sorted.bam -@ 24 -O bam
samtools index ./mapped_reads/$LIBRARY.rm.sorted.bam

###Inspect all results before moving on to variant calling###
