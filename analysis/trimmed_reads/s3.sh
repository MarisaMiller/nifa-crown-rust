#!/bin/bash -l

#PBS -l nodes=1:ppn=1,walltime=36:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N s3_upload

cd $PBS_O_WORKDIR

#This script is to prepare everything for s3 upload and removal from MSI primary storage
tar cvzf trimmed_data.tar.gz *E.gz

s3cmd put -P trimmed_data.tar.gz s3://nifa_crown_rust

rm *.gz
