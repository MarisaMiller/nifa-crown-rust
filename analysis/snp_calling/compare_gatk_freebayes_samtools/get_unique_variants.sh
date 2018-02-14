#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=60:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N unique_variants

cd $PBS_O_WORKDIR

module load liblzma/5.2.2 gcc/4.9.2

#Install your own vcflib for filtering because the version on MSI is old and broken
VCFLIB="/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin"
GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"

#Use decomposed variants as input to get unique variants for each caller
$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 2015_isolates.gatk.filter.sites.vcf 2015_isolates.freebayes.filter_for_compare.sites.vcf > 2015_isolates_unique_to_freebayesVSgatk.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 2015_isolates.samtools.filter.sites.vcf 2015_isolates_unique_to_freebayesVSgatk.vcf > 2015_isolates_unique_to_freebayesVSall.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 2015_isolates.freebayes.filter_for_compare.sites.vcf 2015_isolates.gatk.filter.sites.vcf > 2015_isolates_unique_to_gatkVSfreebayes.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 2015_isolates.samtools.filter.sites.vcf 2015_isolates_unique_to_gatkVSfreebayes.vcf > 2015_isolates_unique_to_gatkVSall.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 2015_isolates.freebayes.filter_for_compare.sites.vcf 2015_isolates.samtools.filter.sites.vcf > 2015_isolates_unique_to_samtoolsVSfreebayes.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 2015_isolates.gatk.filter.sites.vcf 2015_isolates_unique_to_samtoolsVSfreebayes.vcf > 2015_isolates_unique_to_samtoolsVSall.vcf

cat 2015_isolates_unique_to_freebayesVSall.vcf | $VCFLIB/vcfstats -l > 2015_isolate_freebayes_unique_compare_stats.txt
cat 2015_isolates_unique_to_gatkVSall.vcf | $VCFLIB/vcfstats -l > 2015_isolate_gatk_unique_compare_stats.txt
cat 2015_isolates_unique_to_samtoolsVSall.vcf | $VCFLIB/vcfstats -l > 2015_isolate_samtools_unique_compare_stats.txt

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 1990_isolates.gatk.filter.sites.vcf 1990_isolates.freebayes.filter_for_compare.sites.vcf > 1990_isolates_unique_to_freebayesVSgatk.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 1990_isolates.samtools.filter.sites.vcf 1990_isolates_unique_to_freebayesVSgatk.vcf > 1990_isolates_unique_to_freebayesVSall.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 1990_isolates.freebayes.filter_for_compare.sites.vcf 1990_isolates.gatk.filter.sites.vcf > 1990_isolates_unique_to_gatkVSfreebayes.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 1990_isolates.samtools.filter.sites.vcf 1990_isolates_unique_to_gatkVSfreebayes.vcf > 1990_isolates_unique_to_gatkVSall.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 1990_isolates.freebayes.filter_for_compare.sites.vcf 1990_isolates.samtools.filter.sites.vcf > 1990_isolates_unique_to_samtoolsVSfreebayes.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 1990_isolates.gatk.filter.sites.vcf 1990_isolates_unique_to_samtoolsVSfreebayes.vcf > 1990_isolates_unique_to_samtoolsVSall.vcf

cat 1990_isolates_unique_to_freebayesVSall.vcf | $VCFLIB/vcfstats -l > 1990_isolate_freebayes_unique_compare_stats.txt
cat 1990_isolates_unique_to_gatkVSall.vcf | $VCFLIB/vcfstats -l > 1990_isolate_gatk_unique_compare_stats.txt
cat 1990_isolates_unique_to_samtoolsVSall.vcf | $VCFLIB/vcfstats -l > 1990_isolate_samtools_unique_compare_stats.txt
