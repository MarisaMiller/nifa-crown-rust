#!/bin/bash -l

#PBS -l nodes=1:ppn=8,walltime=90:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N calc_R2

cd $PBS_O_WORKDIR

#Note that vcftools in installed and in my $PATH

#Use vcftools to calculate r squared values for 1990 and 2015 years
#--ld-window-bp defines the maximum number of physical bases between the SNPs being tested for LD.
#SNPs were thinned to every 1kb to reduce the comparison size a bit
#Isolates were also split by alternate host region
vcftools --gzvcf ../snp_calling/freebayes/1990_isolates.filter.vcf.gz --thin 1000 --geno-r2 --ld-window-bp 200000 --out 1990_isolates.LD.200kb
vcftools --gzvcf ../snp_calling/freebayes/1990_isolates.filter.vcf.gz --keep mn_bt_list.txt --thin 1000 --geno-r2 --ld-window-bp 200000 --out 1990_isolates.LD.200kb.MNBT
vcftools --gzvcf ../snp_calling/freebayes/1990_isolates.filter.vcf.gz --keep bt_list.txt --thin 1000 --geno-r2 --ld-window-bp 200000 --out 1990_isolates.LD.200kb.BT
vcftools --gzvcf ../snp_calling/freebayes/1990_isolates.filter.vcf.gz --keep no_bt_list.txt --thin 1000 --geno-r2 --ld-window-bp 200000 --out 1990_isolates.LD.200kb.NOBT

vcftools --gzvcf ../snp_calling/freebayes/2015_isolates.filter.ref_removed.vcf.gz --thin 1000 --geno-r2 --ld-window-bp 200000 --out 2015_isolates.LD.200kb
vcftools --gzvcf ../snp_calling/freebayes/2015_isolates.filter.ref_removed.vcf.gz --keep mn_bt_list.txt --thin 1000 --geno-r2 --ld-window-bp 200000 --out 2015_isolates.LD.200kb.MNBT
vcftools --gzvcf ../snp_calling/freebayes/2015_isolates.filter.ref_removed.vcf.gz --keep bt_list.txt --thin 1000 --geno-r2 --ld-window-bp 200000 --out 2015_isolates.LD.200kb.BT
vcftools --gzvcf ../snp_calling/freebayes/2015_isolates.filter.ref_removed.vcf.gz --keep no_bt_list.txt --thin 1000 --geno-r2 --ld-window-bp 200000 --out 2015_isolates.LD.200kb.NOBT
