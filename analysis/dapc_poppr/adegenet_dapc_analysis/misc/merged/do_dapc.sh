#!/bin/bash -l

#PBS -l nodes=1:ppn=15,walltime=4:00:00
#PBS -m abe
#PBS -A kianians
#PBS -M millerme@umn.edu
#PBS -N dapc

cd $PBS_O_WORKDIR

module load R/3.4.0
module load gcc/7.2.0

R CMD BATCH dapc_script.R
