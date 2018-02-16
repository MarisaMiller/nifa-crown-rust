#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=48:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N compare_variants

cd $PBS_O_WORKDIR

module load liblzma/5.2.2 gcc/4.9.2

#Install your own vcflib for filtering because the version on MSI is old and broken
VCFLIB="/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin"
GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"

#Decompose variants to get sites

cat ../gatk/2015_isolates.filter.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 | $VCFLIB/vcfstreamsort > 2015_isolates.gatk.filter.sites.vcf &

cat ../gatk/1990_isolates.filter.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 | $VCFLIB/vcfstreamsort > 1990_isolates.gatk.filter.sites.vcf &

cat ../freebayes/2015_isolates.filter_for_compare.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 | $VCFLIB/vcfstreamsort > 2015_isolates.freebayes.filter_for_compare.sites.vcf &

cat ../freebayes/1990_isolates.filter_for_compare.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 | $VCFLIB/vcfstreamsort > 1990_isolates.freebayes.filter_for_compare.sites.vcf &

cat ../samtools/2015_isolates.filter.mod.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 | $VCFLIB/vcfstreamsort > 2015_isolates.samtools.filter.sites.vcf &

cat ../samtools/1990_isolates.filter.mod.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 | $VCFLIB/vcfstreamsort > 1990_isolates.samtools.filter.sites.vcf &

#These two lines are just to compare 2015 and 1990 freebayes results
cat ../freebayes/2015_isolates.filter.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 > 2015_isolates.freebayes.filter.sites.vcf &

cat ../freebayes/1990_isolates.filter.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | cut -f -8 > 1990_isolates.freebayes.filter.sites.vcf &

wait

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.freebayes.filter_for_compare.sites.vcf 2015_isolates.gatk.filter.sites.vcf > 2015_isolates_freebayes_gatk.intersect.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.samtools.filter.sites.vcf 2015_isolates_freebayes_gatk.intersect.vcf > 2015_gatk_freebayes_samtools_isolates.intersect.vcf

cat 2015_isolates.freebayes.filter_for_compare.sites.vcf | $VCFLIB/vcfstats -l > 2015_isolate_freebayes_for_compare_stats.txt
cat 2015_isolates.gatk.filter.sites.vcf | $VCFLIB/vcfstats -l > 2015_isolate_gatk_stats.txt
cat 2015_isolates.samtools.filter.sites.vcf | $VCFLIB/vcfstats -l > 2015_isolate_samtools_stats.txt
cat 2015_isolates_freebayes_gatk.intersect.vcf | $VCFLIB/vcfstats -l > 2015_isolate_freebayes_gatk_stats.txt
$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.freebayes.filter_for_compare.sites.vcf 2015_isolates.samtools.filter.sites.vcf | $VCFLIB/vcfstats -l > 2015_isolate_freebayes_samtools_stats.txt
$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.gatk.filter.sites.vcf 2015_isolates.samtools.filter.sites.vcf | $VCFLIB/vcfstats -l > 2015_isolate_gatk_samtools_stats.txt
cat 2015_gatk_freebayes_samtools_isolates.intersect.vcf | $VCFLIB/vcfstats -l > 2015_gatk_freebayes_samtools_isolate_intersect_stats.txt

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 1990_isolates.freebayes.filter_for_compare.sites.vcf 1990_isolates.gatk.filter.sites.vcf > 1990_isolates_freebayes_gatk.intersect.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 1990_isolates.samtools.filter.sites.vcf 1990_isolates_freebayes_gatk.intersect.vcf > 1990_gatk_freebayes_samtools_isolates.intersect.vcf

cat 1990_isolates.freebayes.filter_for_compare.sites.vcf | $VCFLIB/vcfstats -l > 1990_isolate_freebayes_for_compare_stats.txt
cat 1990_isolates.gatk.filter.sites.vcf | $VCFLIB/vcfstats -l > 1990_isolate_gatk_stats.txt
cat 1990_isolates.samtools.filter.sites.vcf | $VCFLIB/vcfstats -l > 1990_isolate_samtools_stats.txt
cat 1990_isolates_freebayes_gatk.intersect.vcf | $VCFLIB/vcfstats -l > 1990_isolate_freebayes_gatk_stats.txt
$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 1990_isolates.freebayes.filter_for_compare.sites.vcf 1990_isolates.samtools.filter.sites.vcf | $VCFLIB/vcfstats -l > 1990_isolate_freebayes_samtools_stats.txt
$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 1990_isolates.gatk.filter.sites.vcf 1990_isolates.samtools.filter.sites.vcf | $VCFLIB/vcfstats -l > 1990_isolate_gatk_samtools_stats.txt
cat 1990_gatk_freebayes_samtools_isolates.intersect.vcf | $VCFLIB/vcfstats -l > 1990_gatk_freebayes_samtools_isolate_intersect_stats.txt

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.freebayes.filter.sites.vcf 1990_isolates.freebayes.filter.sites.vcf > 1990_2015_isolates.freebayes.intersect.vcf

cat 2015_isolates.freebayes.filter.sites.vcf | $VCFLIB/vcfstats -l > 2015_isolate_freebayes_stats.txt
cat 1990_isolates.freebayes.filter.sites.vcf | $VCFLIB/vcfstats -l > 1990_isolate_freebayes_stats.txt
cat 1990_2015_isolates.freebayes.intersect.vcf | $VCFLIB/vcfstats -l > 1990_2015_isolate_freebayes_intersect_stats.txt

###This is to get final intersection files for side by side analysis comparison with just FreeBayes results###
cat ../gatk/2015_isolates.filter.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | $VCFLIB/vcfstreamsort > 2015_isolates.gatk.filter.primitive.vcf &

cat ../gatk/1990_isolates.filter.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | $VCFLIB/vcfstreamsort > 1990_isolates.gatk.filter.primitive.vcf &

cat ../freebayes/2015_isolates.filter_for_compare.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | $VCFLIB/vcfstreamsort > 2015_isolates.freebayes.filter_for_compare.primitive.vcf &

cat ../freebayes/1990_isolates.filter_for_compare.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | $VCFLIB/vcfstreamsort > 1990_isolates.freebayes.filter_for_compare.primitive.vcf &

cat ../samtools/2015_isolates.filter.mod.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | $VCFLIB/vcfstreamsort > 2015_isolates.samtools.filter.primitive.vcf &

cat ../samtools/1990_isolates.filter.mod.vcf | $VCFLIB/vcfallelicprimitives -k | $VCFLIB/vcffixup - | $VCFLIB/vcfstreamsort > 1990_isolates.samtools.filter.primitive.vcf &

wait

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.freebayes.filter_for_compare.primitive.vcf 2015_isolates.gatk.filter.primitive.vcf | $VCFLIB/vcffixup - > 2015_isolates_freebayes_gatk.intersect.final.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 2015_isolates.samtools.filter.primitive.vcf 2015_isolates_freebayes_gatk.intersect.final.vcf | $VCFLIB/vcffixup - > 2015_gatk_freebayes_samtools_isolates.intersect.final.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 1990_isolates.freebayes.filter_for_compare.primitive.vcf 1990_isolates.gatk.filter.primitive.vcf | $VCFLIB/vcffixup - > 1990_isolates_freebayes_gatk.intersect.final.vcf

$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -i 1990_isolates.samtools.filter.primitive.vcf 1990_isolates_freebayes_gatk.intersect.final.vcf | $VCFLIB/vcffixup - > 1990_gatk_freebayes_samtools_isolates.intersect.final.vcf

#Compress and index files for future analysis with tabix (note, I have tabix and bgzip installed and in my $PATH as these are not on MSI)
bgzip 1990_gatk_freebayes_samtools_isolates.intersect.final.vcf
tabix -p vcf 1990_gatk_freebayes_samtools_isolates.intersect.final.vcf.gz

bgzip 2015_gatk_freebayes_samtools_isolates.intersect.final.vcf
tabix -p vcf 2015_gatk_freebayes_samtools_isolates.intersect.final.vcf.gz
