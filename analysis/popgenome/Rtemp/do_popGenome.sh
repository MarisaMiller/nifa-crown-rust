#!/bin/bash -l

#PBS -l nodes=1:ppn=15,walltime=16:00:00
#PBS -m abe
#PBS -A kianians
#PBS -M millerme@umn.edu
#PBS -N popgenome

cd $PBS_O_WORKDIR

module load R

R CMD BATCH popgenome_script.R
