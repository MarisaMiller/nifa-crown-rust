#!/bin/bash -l

#PBS -l nodes=1:ppn=30,mem=850gb,walltime=94:00:00
#PBS -m abe
#PBS -M millerme@umn.edu
#PBS -N gatk_pipeline

#submit to a higher memory queue and give as long a walltime as possible

cd $PBS_O_WORKDIR

#Load required modules
module load gatk/3.7.0 java liblzma/5.2.2 gcc/4.7.2

export _JAVA_OPTIONS="-Xmx850g"

#The genome was indexed with samtools and a dictionary was created with Picard with the following commands prior to running:
#module load samtools (v1.5)
#samtools faidx Puccinia_coronata_avenae_12SD80.primary.fa
#module load picard (v2.3.0)
#$PTOOL/picard.jar CreateSequenceDictionary R=Puccinia_coronata_avenae_12SD80.primary.fa O=Puccinia_coronata_avenae_12SD80.primary.dict

GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"

#Make an array of input files
input=(/home/kianians/millerme/nifa-crown-rust/analysis/mapped_reads/*.bam)

for ((i=0;i<=${#input[@]};i++))
do
java -jar /panfs/roc/msisoft/gatk/3.7.0/GenomeAnalysisTK.jar -T HaplotypeCaller -R $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa -I "${input[i]}" -o $(basename "${input[i]%.rm.sorted.bam}.raw.g.vcf") -nct 30 --genotyping_mode DISCOVERY --emitRefConfidence GVCF
done

#Then, run GenotypeGVCFs on 2015 and 1990 samples
input_1990_gvcf=(/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/gatk/90*.g.vcf)

input_2015_gvcf=(/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/gatk/1[25]*.g.vcf)

#Get array in proper format with -V option before each sample
GVCF_1990_IN=()
for s in "${input_1990_gvcf[@]}"
do
GVCF_1990_IN+=(-V $s)
done

GVCF_2015_IN=()
for s in "${input_2015_gvcf[@]}"
do
GVCF_2015_IN+=(-V $s)
done

java -jar /panfs/roc/msisoft/gatk/3.7.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa "${GVCF_1990_IN[@]}" -o 1990_isolates.vcf

java -jar /panfs/roc/msisoft/gatk/3.7.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $GENOME/Puccinia_coronata_avenae_12SD80.primary.fa "${GVCF_2015_IN[@]}" -o 2015_isolates.vcf

#Filter SNPs and get stats

VCFLIB="/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin"

cat 1990_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & AC > 0"  > 1990_isolates.filter.vcf

cat 2015_isolates.vcf | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "QUAL > 20 & AC > 0"  > 2015_isolates.filter.vcf

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_het_stats.txt

(for S in `cat 1990_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 1990_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 1990_isolate_filter_hom_stats.txt

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 1" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_het_stats.txt

(for S in `cat 2015_sample_list.txt`; do echo $S; $VCFLIB/vcfkeepsamples 2015_isolates.filter.vcf $S | $VCFLIB/vcffixup - | $VCFLIB/vcffilter -f "AC > 0 & AC = 2" | $VCFLIB/vcfstats -l; echo -en "\\n"; done;) > 2015_isolate_filter_hom_stats.txt
