#!/bin/bash -l

#PBS -l nodes=1:ppn=16,walltime=16:00:00
#PBS -m abe
#PBS -A kianians
#PBS -M millerme@umn.edu
#PBS -N dapc

cd $PBS_O_WORKDIR

module load R/3.4.4

R CMD BATCH dapc_script.R
