#!/bin/bash -l

set -e
set -u
set -o pipefail

# Run this script locally because it gets killed on MSI (no idea why)

# Set these paths to the locations where you saved the required executables
# Note that the gff3toGenePred can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ as described in the Annovar manual at http://annovar.openbioinformatics.org/en/latest/user-guide/gene/
GFF32GENEPRED="/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/variant_annotation/gff3ToGenePred"
RETRIEVE="/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/variant_annotation/annovar/retrieve_seq_from_fasta.pl"
CONVERT="/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/variant_annotation/annovar/convert2annovar.pl"
ANNOVAR="/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/variant_annotation/annovar/annotate_variation.pl"
GENOME="/home/kianians/millerme/nifa-crown-rust/reference_genomes/"
VCF="/home/kianians/millerme/nifa-crown-rust/analysis/snp_calling/freebayes/"

# mkdir for db
mkdir 12SD80_db

# Convert the GFF3 into 'genePred' format for predicting amino acid sequence impact.
${GFF32GENEPRED} ${GENOME}Puccinia_coronata_avenae_12SD80.p.gff3 12SD80_db/12SD80_p_refGene.txt

# Get the transcript sequences from the GTF, used for translating codons etc.
${RETRIEVE} 12SD80_db/12SD80_p_refGene.txt --format refGene --seqfile ${GENOME}Puccinia_coronata_avenae_12SD80.primary.fa --outfile 12SD80_db/12SD80_p_refGeneMrna.fa

# Convert the VCF into an ANNOVAR input file
# The --allsample option will be used for SNP calling in multiple samples at once
# Also, when --withfreq is set, it will print out the allele frequency of each SNP in the VCF file, based on all samples within the file. Because we are not looking at all samples as a whole, the individual genotypes will not be considered here, so the output file should contain all loci from the input file.
${CONVERT} --format vcf4 --allsample --withfreq --includeinfo --outfile 1990_Annovar_in.txt ${VCF}1990_isolates.filter.vcf.gz
${CONVERT} --format vcf4 --allsample --withfreq --includeinfo --outfile 2015_Annovar_in.txt ${VCF}2015_isolates.filter.vcf.gz

# Then annotate the variants
${ANNOVAR} --geneanno --dbtype refGene --separate --buildver 12SD80_p 1990_Annovar_in.txt 12SD80_db/
${ANNOVAR} --geneanno --dbtype refGene --separate --buildver 12SD80_p 2015_Annovar_in.txt 12SD80_db/

# Here I will just parse out synonymous exonic variants, but the same strategy can be used for any variant annotations
(zgrep "^#" ${VCF}1990_isolates.filter.vcf.gz; grep -w "synonymous SNV" 1990_Annovar_in.txt.exonic_variant_function | cut -f 12-) > 1990_isolates.filter.syn.vcf
(zgrep "^#" ${VCF}2015_isolates.filter.vcf.gz; grep -w "synonymous SNV" 2015_Annovar_in.txt.exonic_variant_function | cut -f 12-) > 2015_isolates.filter.syn.vcf
