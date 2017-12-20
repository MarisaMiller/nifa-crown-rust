#!/bin/bash -l

#PBS -l nodes=1:ppn=24,walltime=96:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N kWIP

cd $PBS_O_WORKDIR

#Load required modules
module load python/2.7.9 gcc/4.8.0 bedtools samtools/1.5

source ../../scripts/khmerEnv/bin/activate

#Generate hashes of reads mapped to 12SD80 for each isolate
mkdir reads
mkdir hashes

R=(/home/kianians/millerme/nifa-crown-rust/analysis/mapped_reads/*.bam)

#samtools view with the -F 4 option saves only mapped reads from the bam files, and then samtools fastq converts this to an interleaved fastq file
for ((i=0;i<=${#R[@]};i++))
do
samtools view -@ 24 -u -F 4 -b "${R[i]}" | samtools fastq -@ 24 - > ./reads/$(basename "${R[i]%.rm.sorted.bam}.fq.gz")
load-into-counting.py -N 1 -k 20 -b -T 24 -f -s tsv -M 100e9 ./hashes/$(basename "${R[i]%.rm.sorted.bam}.ct.gz") ./reads/$(basename "${R[i]%.rm.sorted.bam}.fq.gz")
done

#Distance calculation for 2015
kwip                       \
	-t 24              \
	-k 2015.kern       \
	-d 2015.dist       \
	./hashes/1[25]*.ct.gz

#Distance calculation for 1990
kwip                       \
        -t 24              \
        -k 1990.kern       \
        -d 1990.dist       \
        ./hashes/90*.ct.gz

#Distance calculation for 1990 and 2015 combined
kwip                       \
        -t 24              \
        -k 1990_2015.kern       \
        -d 1990_2015.dist       \
        ./hashes/*.ct.gz

#dist and kern files were summarize with accompanying R script
