#module load R/3.4.4

library(PopGenome)

#Have to run each time if using the ff package
#vcf_1990 <- readData("../1990", format="VCF", gffpath="../GFF")
vcf_1990 <- readData("../1990", format="VCF", gffpath="../GFF_1990")
vcf_2015 <- readData("../2015", format="VCF", gffpath="../GFF")

#all_genes_1990 <- splitting.data(vcf_1990, subsites="gene", whole.data = FALSE)
#all_genes_2015 <- splitting.data(vcf_2015, subsites="gene", whole.data = FALSE)

contig_files_1990 <- dput(as.character(list.files("../FASTA_1990", full.names = TRUE)))
contig_files <- dput(as.character(list.files("../FASTA", full.names = TRUE)))

vcf_1990_syn <- set.synnonsyn(vcf_1990, ref.chr = contig_files_1990)
vcf_2015_syn <- set.synnonsyn(vcf_2015, ref.chr = contig_files)

save(vcf_1990_syn, file="vcf_1990_syn.Rdata")
save(vcf_2015_syn, file="vcf_2015_syn.Rdata")

#This doesn't work for split data because object@region.data@Coding.matrix2 is not generated on split data in whole.data = FALSE. However, if whole.data = TRUE splitting.data doesn't finish properly
#all_genes_1990_syn <- set.synnonsyn(all_genes_1990, ref.chr = contig_files_1990)
#all_genes_2015_syn <- set.synnonsyn(all_genes_2015, ref.chr = contig_files)

#save(all_genes_1990_syn, file="all_genes_1990_syn.Rdata")
#save(all_genes_2015_syn, file="all_genes_2015_syn.Rdata")

#I don't think this module works either. It just seems to stall without progressing
#vcf_1990_4_gamete <- recomb.stats(vcf_1990)
#vcf_2015_4_gamete <- recomb.stats(vcf_2015)

#save(vcf_1990_4_gamete, file="vcf_1990_4_gamete.Rdata")
#save(vcf_2015_4_gamete, file="vcf_2015_4_gamete.Rdata")
