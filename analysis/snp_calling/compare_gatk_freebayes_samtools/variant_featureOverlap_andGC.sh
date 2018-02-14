#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=72:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N overlap_variants

cd $PBS_O_WORKDIR

module load bedtools/2.25.0 liblzma/5.2.2 gcc/4.7.2

GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes"
VCFLIB="/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin"

#First separate SNPs and INDELs
cat 1990_gatk_freebayes_samtools_isolates.intersect.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 1990_gatk_freebayes_samtools_isolates.intersect.snp.vcf

cat 2015_gatk_freebayes_samtools_isolates.intersect.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 2015_gatk_freebayes_samtools_isolates.intersect.snp.vcf

cat 1990_isolates_unique_to_freebayesVSall.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 1990_isolates_unique_to_freebayesVSall.snp.vcf

cat 1990_isolates_unique_to_samtoolsVSall.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 1990_isolates_unique_to_samtoolsVSall.snp.vcf

cat 1990_isolates_unique_to_gatkVSall.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 1990_isolates_unique_to_gatkVSall.snp.vcf

cat 2015_isolates_unique_to_freebayesVSall.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 2015_isolates_unique_to_freebayesVSall.snp.vcf

cat 2015_isolates_unique_to_samtoolsVSall.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 2015_isolates_unique_to_samtoolsVSall.snp.vcf

cat 2015_isolates_unique_to_gatkVSall.vcf | $VCFLIB/vcfsnps | $VCFLIB/vcffixup - > 2015_isolates_unique_to_gatkVSall.snp.vcf

cat 1990_gatk_freebayes_samtools_isolates.intersect.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 1990_gatk_freebayes_samtools_isolates.intersect.indel.vcf

cat 2015_gatk_freebayes_samtools_isolates.intersect.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 2015_gatk_freebayes_samtools_isolates.intersect.indel.vcf

cat 1990_isolates_unique_to_freebayesVSall.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 1990_isolates_unique_to_freebayesVSall.indel.vcf

cat 1990_isolates_unique_to_samtoolsVSall.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 1990_isolates_unique_to_samtoolsVSall.indel.vcf

cat 1990_isolates_unique_to_gatkVSall.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 1990_isolates_unique_to_gatkVSall.indel.vcf

cat 2015_isolates_unique_to_freebayesVSall.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 2015_isolates_unique_to_freebayesVSall.indel.vcf

cat 2015_isolates_unique_to_samtoolsVSall.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 2015_isolates_unique_to_samtoolsVSall.indel.vcf

cat 2015_isolates_unique_to_gatkVSall.vcf | $VCFLIB/vcfindels | $VCFLIB/vcffixup - > 2015_isolates_unique_to_gatkVSall.indel.vcf

#Annotate genes and repeats with counts of overlapping variants
bioawk -c gff '{if ($feature == "gene") print $0}' $GENOME/Puccinia_coronata_avenae_12SD80.p.gff3 | sort -k1,1 -k4,4n | bedtools merge -i stdin > $GENOME/Puccinia_coronata_avenae_12SD80.p.genes.merged.bed

cat $GENOME/Puccinia_coronata_avenae_12SD80.repeats.p.gff3 | sort -k1,1 -k4,4n | bedtools merge -i stdin > $GENOME/Puccinia_coronata_avenae_12SD80.p.repeats.merged.bed

bedtools annotate -counts -i $GENOME/Puccinia_coronata_avenae_12SD80.p.genes.merged.bed -files 1990_gatk_freebayes_samtools_isolates.intersect.snp.vcf 2015_gatk_freebayes_samtools_isolates.intersect.snp.vcf 1990_isolates_unique_to_freebayesVSall.snp.vcf 1990_isolates_unique_to_samtoolsVSall.snp.vcf 1990_isolates_unique_to_gatkVSall.snp.vcf 2015_isolates_unique_to_freebayesVSall.snp.vcf 2015_isolates_unique_to_samtoolsVSall.snp.vcf 2015_isolates_unique_to_gatkVSall.snp.vcf -names snp_shared_1990 snp_shared_2015 snp_unique_1990_freebayes snp_unique_1990_samtools snp_unique_1990_gatk snp_unique_2015_freebayes snp_unique_2015_samtools snp_unique_2015_gatk > snp_variants_geneAnnotate.txt

bedtools annotate -counts -i $GENOME/Puccinia_coronata_avenae_12SD80.p.repeats.merged.bed -files 1990_gatk_freebayes_samtools_isolates.intersect.snp.vcf 2015_gatk_freebayes_samtools_isolates.intersect.snp.vcf 1990_isolates_unique_to_freebayesVSall.snp.vcf 1990_isolates_unique_to_samtoolsVSall.snp.vcf 1990_isolates_unique_to_gatkVSall.snp.vcf 2015_isolates_unique_to_freebayesVSall.snp.vcf 2015_isolates_unique_to_samtoolsVSall.snp.vcf 2015_isolates_unique_to_gatkVSall.snp.vcf -names snp_shared_1990 snp_shared_2015 snp_unique_1990_freebayes snp_unique_1990_samtools snp_unique_1990_gatk snp_unique_2015_freebayes snp_unique_2015_samtools snp_unique_2015_gatk > snp_variants_repeatAnnotate.txt

bedtools annotate -counts -i $GENOME/Puccinia_coronata_avenae_12SD80.p.genes.merged.bed -files 1990_gatk_freebayes_samtools_isolates.intersect.indel.vcf 2015_gatk_freebayes_samtools_isolates.intersect.indel.vcf 1990_isolates_unique_to_freebayesVSall.indel.vcf 1990_isolates_unique_to_samtoolsVSall.indel.vcf 1990_isolates_unique_to_gatkVSall.indel.vcf 2015_isolates_unique_to_freebayesVSall.indel.vcf 2015_isolates_unique_to_samtoolsVSall.indel.vcf 2015_isolates_unique_to_gatkVSall.indel.vcf -names indel_shared_1990 indel_shared_2015 indel_unique_1990_freebayes indel_unique_1990_samtools indel_unique_1990_gatk indel_unique_2015_freebayes indel_unique_2015_samtools indel_unique_2015_gatk > indel_variants_geneAnnotate.txt

bedtools annotate -counts -i $GENOME/Puccinia_coronata_avenae_12SD80.p.repeats.merged.bed -files 1990_gatk_freebayes_samtools_isolates.intersect.indel.vcf 2015_gatk_freebayes_samtools_isolates.intersect.indel.vcf 1990_isolates_unique_to_freebayesVSall.indel.vcf 1990_isolates_unique_to_samtoolsVSall.indel.vcf 1990_isolates_unique_to_gatkVSall.indel.vcf 2015_isolates_unique_to_freebayesVSall.indel.vcf 2015_isolates_unique_to_samtoolsVSall.indel.vcf 2015_isolates_unique_to_gatkVSall.indel.vcf -names indel_shared_1990 indel_shared_2015 indel_unique_1990_freebayes indel_unique_1990_samtools indel_unique_1990_gatk indel_unique_2015_freebayes indel_unique_2015_samtools indel_unique_2015_gatk > indel_variants_repeatAnnotate.txt

#Calculate GC content in windows surrounding unique and shared variants
#To calculate correctly, the header lines need to be removed, and then GC content is calculated 1,000bp upstream and downstream of variant start site
grep -v "^#" 1990_gatk_freebayes_samtools_isolates.intersect.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_gatk_freebayes_samtools_isolates.intersect.gc.snp.txt

grep -v "^#" 1990_gatk_freebayes_samtools_isolates.intersect.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_gatk_freebayes_samtools_isolates.intersect.gc.indel.txt

grep -v "^#" 2015_gatk_freebayes_samtools_isolates.intersect.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_gatk_freebayes_samtools_isolates.intersect.gc.snp.txt

grep -v "^#" 2015_gatk_freebayes_samtools_isolates.intersect.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_gatk_freebayes_samtools_isolates.intersect.gc.indel.txt

grep -v "^#" 1990_isolates_unique_to_freebayesVSall.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_isolates_unique_to_freebayesVSall.gc.snp.txt

grep -v "^#" 1990_isolates_unique_to_freebayesVSall.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_isolates_unique_to_freebayesVSall.gc.indel.txt

grep -v "^#" 1990_isolates_unique_to_samtoolsVSall.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_isolates_unique_to_samtoolsVSall.gc.snp.txt

grep -v "^#" 1990_isolates_unique_to_samtoolsVSall.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_isolates_unique_to_samtoolsVSall.gc.indel.txt

grep -v "^#" 1990_isolates_unique_to_gatkVSall.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_isolates_unique_to_gatkVSall.gc.snp.txt

grep -v "^#" 1990_isolates_unique_to_gatkVSall.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 1990_isolates_unique_to_gatkVSall.gc.indel.txt

grep -v "^#" 2015_isolates_unique_to_freebayesVSall.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_isolates_unique_to_freebayesVSall.gc.snp.txt

grep -v "^#" 2015_isolates_unique_to_freebayesVSall.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_isolates_unique_to_freebayesVSall.gc.indel.txt

grep -v "^#" 2015_isolates_unique_to_samtoolsVSall.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_isolates_unique_to_samtoolsVSall.gc.snp.txt

grep -v "^#" 2015_isolates_unique_to_samtoolsVSall.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_isolates_unique_to_samtoolsVSall.gc.indel.txt

grep -v "^#" 2015_isolates_unique_to_gatkVSall.snp.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_isolates_unique_to_gatkVSall.gc.snp.txt

grep -v "^#" 2015_isolates_unique_to_gatkVSall.indel.vcf | awk -v OFS="\\t" '{print $1,$2,$2}' | bedtools slop -i stdin -g $GENOME/Puccinia_coronata_avenae_12SD80.primary.chromlen -b 1000 | bedtools nuc -fi $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -bed stdin > 2015_isolates_unique_to_gatkVSall.gc.indel.txt

#Alternative method to get intersections, although output is larger and harder to work with. 
#Is useful if actualy variant information is needed after intersection
#Intersect unique and shared vcf files with protein_coding genes and repeats
bioawk -c gff '{if ($feature == "gene") print $0}' $GENOME/Puccinia_coronata_avenae_12SD80.p.gff3 | sort -k1,1 -k4,4n | bedtools intersect -wa -wb -a stdin -b 1990_gatk_freebayes_samtools_isolates.intersect.vcf 2015_gatk_freebayes_samtools_isolates.intersect.vcf 1990_isolates_unique_to_freebayesVSall.vcf 1990_isolates_unique_to_samtoolsVSall.vcf 1990_isolates_unique_to_gatkVSall.vcf 2015_isolates_unique_to_freebayesVSall.vcf 2015_isolates_unique_to_samtoolsVSall.vcf 2015_isolates_unique_to_gatkVSall.vcf -names shared_1990 shared_2015 unique_1990_freebayes unique_1990_samtools unique_1990_gatk unique_2015_freebayes unique_2015_samtools unique_2015_gatk > variants_geneIntersect.txt

bedtools intersect -wa -wb -a $GENOME/Puccinia_coronata_avenae_12SD80.repeats.p.gff3 -b 1990_gatk_freebayes_samtools_isolates.intersect.vcf 2015_gatk_freebayes_samtools_isolates.intersect.vcf 1990_isolates_unique_to_freebayesVSall.vcf 1990_isolates_unique_to_samtoolsVSall.vcf 1990_isolates_unique_to_gatkVSall.vcf 2015_isolates_unique_to_freebayesVSall.vcf 2015_isolates_unique_to_samtoolsVSall.vcf 2015_isolates_unique_to_gatkVSall.vcf -names shared_1990 shared_2015 unique_1990_freebayes unique_1990_samtools unique_1990_gatk unique_2015_freebayes unique_2015_samtools unique_2015_gatk > variants_repeatIntersect.txt
