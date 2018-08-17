#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=4:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N variants_effectors

cd $PBS_O_WORKDIR

module load bedtools/2.25.0 #2.27 produces an error with this intersect command. Do not use!

#Identify annotated variants overlapping effectors

#The first awk command is a strategy to "trick" bedtools into overlapping a file (Annovar output) that's not quite a vcf file with a bed file. The command moves the first 11 columns to the end of the vcf, and since bedtools only need the coordinates from the vcf file to be in the right place, the rest of the file gets dragged along in the intersection. I used the strategy at this link (https://stackoverflow.com/questions/26892639/move-certain-columns-to-end-using-awk). However, I couldn't get the first column not to print, so I piped the output into cut to remove the lingering first column and then into sed to remove trailing tab. I also keep the vcf header with the grep command.

#Also, when this job is submitted on MSI, make sure to escape the slashes in front of \t with awk (so it's \\t). If run on the command line with plain bash shell, this is not needed.

(zgrep "^#" ../../snp_calling/freebayes/1990_isolates.filter.vcf.gz; awk -F"\\t" -v a=1 -v b=12 '{for (i=1;i<=NF;i+=i==a?b-a:1) {printf "%s\\t",$i};for (i=a;i<b;i++) {printf "%s\\t",$i};print""}' ../../snp_calling/variant_annotation/1990_Annovar_in.txt.exonic_variant_function | cut -f 2- | sed 's/\t$//') | bedtools intersect -wa -wb -a stdin -b ../pav_analysis/effector_genes_sd80_p.sorted.gff3 > 1990_annovar_exonic_effector_variants.txt

(zgrep "^#" ../../snp_calling/freebayes/2015_isolates.filter.vcf; awk -F"\\t" -v a=1 -v b=12 '{for (i=1;i<=NF;i+=i==a?b-a:1) {printf "%s\\t",$i};for (i=a;i<b;i++) {printf "%s\\t",$i};print""}' ../../snp_calling/variant_annotation/2015_Annovar_in.txt.exonic_variant_function | cut -f 2- | sed 's/\t$//') | bedtools intersect -wa -wb -a stdin -b ../pav_analysis/effector_genes_sd80_p.sorted.gff3 > 2015_annovar_exonic_effector_variants.txt
