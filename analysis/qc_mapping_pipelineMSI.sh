#!/bin/bash -l

#PBS -l nodes=1:ppn=24,walltime=24:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N SNP_pipeline

cd $PBS_O_WORKDIR

#Load required modules
module load trimmomatic fastx_toolkit bwa/0.7.15 samtools/1.5

RAW_DATA="/home/kianians/millerme/nifa-crown-rust/raw_data"

#Make array of filenames
filelist=($RAW_DATA/*.fastq)

#Generate an array with name after _ gone
filelist2=()
for F in "${filelist[@]}"
do
filelist2+=("${F%_*}")
done

#Remove duplicate values in list to just get unique basenames
uniq=($(printf "%s\n" "${filelist2[@]}" | uniq)); echo "${uniq[@]}"

#Then, use the unique array above in shell script w/ SNP pipeline to submit to MSI queues
for F in "${uniq[@]}"
do
bash ../scripts/QC_mapping_pipeline.sh -n ../raw_reads/$F -g ../reference_genomes/Puccinia_coronata_avenae_12SD80.primary.fa
done
