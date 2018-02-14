#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=96:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N samtools

cd $PBS_O_WORKDIR

#Load required modules
module load samtools liblzma/5.2.2 gcc/4.7.2 bcftools/1.6

#Path to genome file locations. Make sure the reference is indexed with samtools faidx and it is in the same directory.
GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"

#Install your own vcflib for filtering because the version on MSI is old and broken
VCFLIB="/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin"

#Call variants with mpileup
samtools mpileup -b 1990_bam_list.txt -f $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -u -v | bcftools call --multiallelic-caller --variants-only -Ov -o 1990_isolates.vcf &

samtools mpileup -b 2015_bam_list.txt -f $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -u -v | bcftools call --multiallelic-caller --variants-only -Ov -o 2015_isolates.vcf &

wait

#Filter variants and get stats

cat 1990_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & AC > 0"  > 1990_isolates.filter.vcf &

cat 2015_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & AC > 0"  > 2015_isolates.filter.vcf &

wait

#Fix lowercase bases in samtools VCF records (causes issues with future comparison steps)
(grep "^#" 1990_isolates.filter.vcf; grep -v "^#" 1990_isolates.filter.vcf | awk '{printf "%s\\t%s\\t%s\\t%s\\t%s\\t",$1,$2,$3,toupper($4),toupper($5); for(i=6; i<=NF; i++){printf("%s%s", $(i), i<NF ? "\\t" : "\\n")}}') > 1990_isolates.filter.mod.vcf

(grep "^#" 2015_isolates.filter.vcf; grep -v "^#" 2015_isolates.filter.vcf | awk '{printf "%s\\t%s\\t%s\\t%s\\t%s\\t",$1,$2,$3,toupper($4),toupper($5); for(i=6; i<=NF; i++){printf("%s%s", $(i), i<NF ? "\\t" : "\\n")}}') > 2015_isolates.filter.mod.vcf

#Make a sample list with sample names exactly as they appear in the vcf file separated by a single space
#Then, fixup the alternate allele counts for each individual sample, filter out sites that are the same as the ref for a given isolate, and then count het and hom sites
(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.mod.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_het_stats.txt

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.mod.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_hom_stats.txt

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.mod.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_het_stats.txt
 
(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.mod.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_hom_stats.txt

#For filtering out only SNPs
#cat 1990_isolates.filter.vcf | $VCFLIB/vcffilter -f "! ( INDEL )" | $VCFLIB/vcffixup - > 1990_isolates.snps.filter.vcf

#cat 2015_isolates.filter.vcf | $VCFLIB/vcffilter -f "! ( INDEL )" | $VCFLIB/vcffixup - > 2015_isolates.snps.filter.vcf

#(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.snps.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_snps_filter_het_stats.txt

#(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.snps.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_snps_filter_hom_stats.txt

#(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.snps.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_snps_filter_het_stats.txt
 
#(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.snps.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_snps_filter_hom_stats.txt
