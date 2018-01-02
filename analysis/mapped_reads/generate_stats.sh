#!/bin/bash -l

#PBS -l nodes=1:ppn=24,walltime=24:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N stats

module load samtools

cd $PBS_O_WORKDIR

(for F in `ls 15*.bam`; do echo $F; samtools flagstat -@24 $F; echo -en "\\n"; done;) > mapping_stats_2015.txt

(for F in `ls 90*.bam`; do echo $F; samtools flagstat -@24 $F; echo -en "\\n"; done;) > mapping_stats_1990.txt
