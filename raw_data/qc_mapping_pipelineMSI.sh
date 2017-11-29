#!/bin/bash -l

#PBS -l nodes=1:ppn=16,walltime=56:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N qc_mapping_pipeline

#Run this script from inside of the raw_data directory
cd $PBS_O_WORKDIR

#Load required modules
module load trimmomatic fastx_toolkit bwa/0.7.15 samtools/1.5 gnuplot

ANALYSIS="/home/kianians/millerme/nifa-crown-rust/analysis"

#Make output dir for trimmed reads
mkdir $ANALYSIS/trimmed_reads 

#Make output dirs for results of fastx toolkit
mkdir -p $ANALYSIS/trimmed_reads/{statistics,plots/quality_plots,plots/nt_distribution}

#Output directory for mapped reads
mkdir $ANALYSIS/mapped_reads

#Make array of filenames
filelist=(*.fastq.gz)

#Generate an array with name after _ gone
filelist2=()
for F in "${filelist[@]}"
do
filelist2+=("${F%_*}")
done

#Remove duplicate values in list to just get unique basenames
uniq=($(printf "%s\\n" "${filelist2[@]}" | uniq)); echo "${uniq[@]}"

#Then use the unique array above in shell script w/ SNP pipeline to submit to MSI queues

for F in "${uniq[@]}"
do
bash ../scripts/QC_mapping_pipeline.sh -n $F -g ../reference_genomes/Puccinia_coronata_avenae_12SD80.primary.fa -t 16
done
