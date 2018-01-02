#!/bin/bash -l
#PBS -l walltime=36:00:00,nodes=1:ppn=24
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N coverage

cd $PBS_O_WORKDIR

module load bedtools #v2.25
#This script also requires bioawk. I installed this myself on MSI.

#Generate a genome file and a bed file for the reference
GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"
cat $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa | bioawk -c fastx '{ print $name, length($seq) }' |  sort -k1,1 > $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen
awk '{print $1"\t" 0 "\t" $2}' $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen > 12SD80_p.bed

#Make an array of input files
input=(/home/kianians/millerme/nifa-crown-rust/analysis/mapped_reads/*.bam)

#Convert bam files to bed format and then sort.
#Then, use bedtools coverage to calculate the average coverage across each contig in the genome.
for ((i=0;i<=${#input[@]};i++))
do
bedtools bamtobed -i "${input[i]}" | sort -k1,1 -k2,2n -S 75% | bedtools coverage -mean -sorted -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -a 12SD80_p.bed -b - > $(basename "${input[i]%.rm.sorted.bam}.cov")
done
