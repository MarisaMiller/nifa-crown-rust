#Usage: bash QC_mapping_pipeline.sh -n basename_of_each_library_fastq_file -g path/to/genome.fasta -t threads
#The pipeline doesn't need to be modified if your basenames fit the format of XXXX_R1.fastq or XXXX_R2.fastq (where the XXXX can be anything as long as it doesn't contain spaces or underscores, the only underscore should be before the R1 or R2).
#If you have a different file format, then you need to modify the MSI job submission script (placed in the analysis directory).
#Use the rename.py script in the scripts directory to rename files from UMGC to fit this format.
#Make sure your referenced genome in the reference_genomes folder is indexed with BWA and samtools. Here are the commands: bwa index reference.fasta and samtools faidx reference.fasta
#Read through script carefully before running to ensure all assumptions and settings match your system.

#Make a project directory (named whatever you like), and then set up the following subdirectory structures for this pipeline:
#./raw_data    Place the MSI job script and submit it from inside this directory. Raw data can be compressed (modify script below if not).
#./reference_genomes
#./scripts     Add this directory to your $PATH if desired
#./analysis
 
#Check that the correct parameters were given, and if so assign the variables of LIBRARY and GENOME and THREADS
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
		-t|--threads) #Give the number of threads
                THREADS="$2"
                shift
                ;;
		#Exit here if -n, -g, or -t options were not given correctly
		*) echo "Usage: QC_mapping_pipeline.sh -n path/to/filename_prefix -g path/to/genome.fasta -t threads. Unknown option: $1 " >&2 ; exit 1
		;;
	esac
	shift
done

#Give error and exit if the variables were not set correctly
if [ -z ${LIBRARY} ]; then echo "No library directory given! Usage: QC_mapping_pipeline.sh -n path/to/filename_prefix -g path/to/genome.fasta -t threads"; exit 1; fi
if [ -z $GENOME ]; then echo "No reference genome given! Usage: QC_mapping_pipeline.sh -n path/to/filename_prefix -g path/to/genome.fasta -t threads"; exit 1 ; fi

#Print the values assigned to the variables
echo "LIBRARY" = "${LIBRARY}"
echo "GENOME" = "${GENOME}"
echo "THREADS" = "${THREADS}"

#Identify the raw read files. 
fasR=$(ls ${LIBRARY}_R1.fastq.gz)
fasL=$(ls ${LIBRARY}_R2.fastq.gz)

#If input files are not found then exit.
if [ ! -f $fasR ] || [ ! -f $fasL ]
then
	echo "Reads from "${LIBRARY}" could not be found"; exit 1
fi


####################################
######       SNP calling      ###### 
####################################


#Perform filtering and trimming of the raw reads, followed by QC
                
#$TRIMMOMATIC is an MSI specific environment variable. On your system, just replace with the path to the Trimmomatic executables.
#Change the threads to match your system
#Change the adapter and other parameters to fit your requirements

java -jar $TRIMMOMATIC/trimmomatic.jar PE -threads $THREADS $fasR $fasL ../analysis/trimmed_reads/${LIBRARY}_trim_R1.PE ../analysis/trimmed_reads/${LIBRARY}_trim_R1.SE ../analysis/trimmed_reads/${LIBRARY}_trim_R2.PE ../analysis/trimmed_reads/${LIBRARY}_trim_R2.SE ILLUMINACLIP:$TRIMMOMATIC/adapters/all_illumina_adapters.fa:2:30:10

#There should be very few SE reads since our trimming was not very aggressive.
#Whatever SE files are not used into mapping.
fastx_quality_stats -i ../analysis/trimmed_reads/${LIBRARY}_trim_R1.PE -o ../analysis/trimmed_reads/statistics/${LIBRARY}_R1_PE_stats.txt
fastx_quality_stats -i ../analysis/trimmed_reads/${LIBRARY}_trim_R1.SE -o ../analysis/trimmed_reads/statistics/${LIBRARY}_R1_SE_stats.txt
fastx_quality_stats -i ../analysis/trimmed_reads/${LIBRARY}_trim_R2.PE -o ../analysis/trimmed_reads/statistics/${LIBRARY}_R2_PE_stats.txt
fastx_quality_stats -i ../analysis/trimmed_reads/${LIBRARY}_trim_R2.SE -o ../analysis/trimmed_reads/statistics/${LIBRARY}_R2_SE_stats.txt

fastq_quality_boxplot_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R1_PE_stats.txt -o ../analysis/trimmed_reads/plots/quality_plots/${LIBRARY}_R1_PE_quality.png
fastq_quality_boxplot_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R1_SE_stats.txt -o ../analysis/trimmed_reads/plots/quality_plots/${LIBRARY}_R1_SE_quality.png
fastq_quality_boxplot_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R2_PE_stats.txt -o ../analysis/trimmed_reads/plots/quality_plots/${LIBRARY}_R2_PE_quality.png
fastq_quality_boxplot_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R2_SE_stats.txt -o ../analysis/trimmed_reads/plots/quality_plots/${LIBRARY}_R2_SE_quality.png
                
fastx_nucleotide_distribution_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R1_PE_stats.txt -o ../analysis/trimmed_reads/plots/nt_distribution/${LIBRARY}_R1_PE_nt_distr.png
fastx_nucleotide_distribution_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R1_SE_stats.txt -o ../analysis/trimmed_reads/plots/nt_distribution/${LIBRARY}_R1_SE_nt_distr.png
fastx_nucleotide_distribution_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R2_PE_stats.txt -o ../analysis/trimmed_reads/plots/nt_distribution/${LIBRARY}_R2_PE_nt_distr.png
fastx_nucleotide_distribution_graph.sh -i ../analysis/trimmed_reads/statistics/${LIBRARY}_R2_SE_stats.txt -o ../analysis/trimmed_reads/plots/nt_distribution/${LIBRARY}_R2_SE_nt_distr.png

#Align against the reference genome (previously indexed), and then convert to bam, remove duplicates, sort and index.

#Assign unique readgroups with the -R flag so that samples can be used in pooled SNP calling later
bwa mem -t $THREADS -R "@RG\tID:${LIBRARY}\tSM:${LIBRARY}" $GENOME ../analysis/trimmed_reads/${LIBRARY}_trim_R1.PE ../analysis/trimmed_reads/${LIBRARY}_trim_R2.PE | samtools view -bh -@ $THREADS - | samtools sort -@ $THREADS -O bam -T temp - | samtools rmdup - ../analysis/mapped_reads/${LIBRARY}.rm.sorted.bam
samtools index ../analysis/mapped_reads/${LIBRARY}.rm.sorted.bam

###Inspect all results before moving on to variant calling###
