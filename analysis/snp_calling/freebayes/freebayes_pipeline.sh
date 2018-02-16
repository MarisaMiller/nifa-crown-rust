#!/bin/bash -l

#PBS -l nodes=1:ppn=4,walltime=96:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N freebayes

cd $PBS_O_WORKDIR

#Load required modules
module load freebayes/20161103 parallel liblzma/5.2.2 gcc/4.7.2

#Install your own vcflib for filtering because the version on MSI is old and broken
VCFLIB="/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin"

GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"

#Call SNPs with FreeBayes

#Just need to create regions once, comment this out after performing once for a given genome
#/soft/freebayes/20161103/bin/fasta_generate_regions.py $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa 100000 > $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa.regions

#Make an array of the bam alignment files for 2015 and 2012 isolates
input_2015=$(ls /home/kianians/millerme/nifa-crown-rust/analysis/mapped_reads/1[25]*.bam)

#Or for just 1990 samples
input_1990=$(ls /home/kianians/millerme/nifa-crown-rust/analysis/mapped_reads/90*.bam)

#Call SNPs for all 2015 or 1990 isolates at once, and then get stats.

./freebayes-parallel $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa.regions 24 -f $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa $input_2015 > 2015_isolates.vcf

./freebayes-parallel $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa.regions 24 -f $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa $input_1990 > 1990_isolates.vcf

#Filter SNPs and get stats

cat 2015_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & AC > 0"  > 2015_isolates.filter.vcf

cat 1990_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & AC > 0"  > 1990_isolates.filter.vcf

#Make a sample list with sample names exactly as they appear in the vcf file separated by a single space
#Then, fixup the alternate allele counts for each individual sample, filter out sites that are the same as the ref for a given isolate, and then count het and hom sites
(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_stats.txt

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_het_stats.txt

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_hom_stats.txt

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_stats.txt

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_het_stats.txt
 
(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_hom_stats.txt

#This part is for generating stats for comparison with GATK
cat 2015_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & AC > 0"  > 2015_isolates.filter_for_compare.vcf
cat 1990_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & AC > 0"  > 1990_isolates.filter_for_compare.vcf

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter_for_compare.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_for_compare_het_stats.txt

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter_for_compare.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_for_compare_hom_stats.txt

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter_for_compare.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_for_compare_het_stats.txt

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter_for_compare.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_for_compare_hom_stats.txt

#Combine 1990 and 2015 sets together for later analysis with DAPC
$VCFLIB/vcfcombine 2015_isolates.filter.vcf 1990_isolates.filter.vcf > 1990_2015_isolates.filter.vcf

#Compress and index files for future analysis with tabix (note, I have tabix and bgzip installed and in my $PATH as these are not on MSI)
bgzip 1990_isolates.filter.vcf
tabix -p vcf 1990_isolates.filter.vcf.gz

bgzip 2015_isolates.filter.vcf
tabix -p vcf 2015_isolates.filter.vcf.gz

###After thinking more about the issue of filtering out sites that are also het in the reference, I decided that by throwing away variant sites that are the same (het) as the reference, we are actually losing a lot of variation among the samples. Really, we don't care that a given variant is the same or different as 12SD80, but how the samples differ amongst each other with respect to a given position.
###These commands were kept as a log of what I tried with the 2015 isolates as a test.

#Then, do filtering against 12SD80 and get stats

#First, filter out just 12SD80, and do a little further filtering to remove not heterozygous (aka heterokaryotic) sites (as we don't expect those when mapping an isolate to itself)
#$VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf 12SD80 | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AB > 0.2 & AB < 0.8" > 12SD80_filter.vcf

######This gave the same result as above
######$VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf 12SD80 | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC < 2" > 12SD80_filter_alt.vcf

#$VCFLIB/vcfstats -l 12SD80_filter.vcf > 12SD80_filter_stats.txt
#$VCFLIB/vcffilter -f "AC = 1" 12SD80_filter.vcf | $VCFLIB/vcfstats -l > 12SD80_filter_het_stats.txt
#$VCFLIB/vcffilter -f "AC = 2" 12SD80_filter.vcf | $VCFLIB/vcfstats -l > 12SD80_filter_hom_stats.txt

#######Also, removing 12SD80, and then getting stats, did not influence the variant statistic results.
#$VCFLIB/vcfremovesamples 2015_isolates.filter.vcf 12SD80 | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0" > 2015_isolates_no12SD80.filter.vcf

#(for S in `cat sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates_no12SD80.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_no12SD80_filter_het_stats.txt

#(for S in `cat sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates_no12SD80.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_no12SD80_filter_hom_stats.txt

#$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 12SD80_filter.vcf 2015_isolates_no12SD80.filter.vcf > 2015_isolates.ref_filter.vcf

######$VCFLIB/vcfintersect -r $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -v -i 12SD80_filter.vcf 2015_isolates.filter.vcf > 2015_isolates.ref_filter.vcf

#Stats of final 12SD80 filtered multisample.vcf
#(for S in `cat sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.ref_filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_ref_filter_stats.txt

#(for S in `cat sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.ref_filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_ref_filter_het_stats.txt

#(for S in `cat sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.ref_filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_ref_filter_hom_stats.txt
