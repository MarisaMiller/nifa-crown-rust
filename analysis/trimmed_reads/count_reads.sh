#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=24:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N count_reads

module load parallel

cd $PBS_O_WORKDIR

parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: *R1.PE.gz
