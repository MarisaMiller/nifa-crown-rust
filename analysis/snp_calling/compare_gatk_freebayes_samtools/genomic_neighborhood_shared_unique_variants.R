#In this script, the genomic neighborhood (genes vs repeats) of shared and unique variants between different callers is explored

library(dplyr)
library(purrr)
library(ggplot2)

#Set working directory to appropriate location
setwd("~/Documents/Work/Log Book Data/2016/2016_crownrust_NIFA_project/sequencing_data_analysis/gatk_freebayes_samtools_compare/neighborhood_analysis")

#Read in datasets with variant counts in genes and repeats
snp_anno_gene_df <- read.delim("snp_variants_geneAnnotate.txt", header=TRUE, col.names = c("contig", "start", "end", "Shared SNPs 1990", "Shared SNPs 2015", "Unique SNPs 1990 Freebayes", "Unique SNPs 1990 Samtools", "Unique SNPs 1990 GATK", "Unique SNPs 2015 Freebayes", "Unique SNPs 2015 Samtools", "Unique SNPs 2015 GATK"))

indel_anno_gene_df <- read.delim("indel_variants_geneAnnotate.txt", header=TRUE, col.names = c("contig", "start", "end", "Shared INDELs 1990", "Shared INDELs 2015", "Unique INDELs 1990 Freebayes", "Unique INDELs 1990 Samtools", "Unique INDELs 1990 GATK", "Unique INDELs 2015 Freebayes", "Unique INDELs 2015 Samtools", "Unique INDELs 2015 GATK"))

snp_anno_repeat_df <- read.delim("snp_variants_repeatAnnotate.txt", header=TRUE, col.names = c("contig", "start", "end", "Shared SNPs 1990", "Shared SNPs 2015", "Unique SNPs 1990 Freebayes", "Unique SNPs 1990 Samtools", "Unique SNPs 1990 GATK", "Unique SNPs 2015 Freebayes", "Unique SNPs 2015 Samtools", "Unique SNPs 2015 GATK"))

indel_anno_repeat_df <- read.delim("indel_variants_repeatAnnotate.txt", header=TRUE, col.names = c("contig", "start", "end", "Shared INDELs 1990", "Shared INDELs 2015", "Unique INDELs 1990 Freebayes", "Unique INDELs 1990 Samtools", "Unique INDELs 1990 GATK", "Unique INDELs 2015 Freebayes", "Unique INDELs 2015 Samtools", "Unique INDELs 2015 GATK"))

#Read in the summary of total unique and shared variants
summary_df <- read.delim("../gatk_freebayes_samtools_compare_1990_2015.txt", header=TRUE)

#Pull out the relevant values for normalizing the variant counts
total_int_1990_snp <- summary_df[11,4]
total_int_2015_snp <- summary_df[4,4]
freebayes_1990_unique_snps <- summary_df[9,4] - summary_df[11,4] - (summary_df[12,4] - summary_df[11,4]) - (summary_df[14,4] -  summary_df[11,4])
samtools_1990_unique_snps <- summary_df[10,4] - summary_df[11,4] - (summary_df[13,4] - summary_df[11,4]) - (summary_df[14,4] -  summary_df[11,4])
gatk_1990_unique_snps <- summary_df[8,4] - summary_df[11,4] - (summary_df[12,4] - summary_df[11,4]) - (summary_df[13,4] -  summary_df[11,4])
freebayes_2015_unique_snps <- summary_df[2,4] - summary_df[4,4] - (summary_df[5,4] - summary_df[4,4]) - (summary_df[7,4] -  summary_df[4,4])
samtools_2015_unique_snps <- summary_df[3,4] - summary_df[4,4] - (summary_df[6,4] - summary_df[4,4]) - (summary_df[7,4] -  summary_df[4,4])
gatk_2015_unique_snps <- summary_df[1,4] - summary_df[4,4] - (summary_df[5,4] - summary_df[4,4]) - (summary_df[6,4] -  summary_df[4,4])
total_int_1990_indel <- summary_df[11,5]
total_int_2015_indel <- summary_df[4,5]
freebayes_1990_unique_indels <- summary_df[9,5] - summary_df[11,5] - (summary_df[12,5] - summary_df[11,5]) - (summary_df[14,5] -  summary_df[11,5])
samtools_1990_unique_indels <- summary_df[10,5] - summary_df[11,5] - (summary_df[13,5] - summary_df[11,5]) - (summary_df[14,5] -  summary_df[11,5])
gatk_1990_unique_indels <- summary_df[8,5] - summary_df[11,5] - (summary_df[12,5] - summary_df[11,5]) - (summary_df[13,5] -  summary_df[11,5])
freebayes_2015_unique_indels <- summary_df[2,5] - summary_df[4,5] - (summary_df[5,5] - summary_df[4,5]) - (summary_df[7,5] -  summary_df[4,5])
samtools_2015_unique_indels <- summary_df[3,5] - summary_df[4,5] - (summary_df[6,5] - summary_df[4,5]) - (summary_df[7,5] -  summary_df[4,5])
gatk_2015_unique_indels <- summary_df[1,5] - summary_df[4,5] - (summary_df[5,5] - summary_df[4,5]) - (summary_df[6,5] -  summary_df[4,5])

#Gather values for normalizing
norm_snps <- c(total_int_1990_snp, total_int_2015_snp, freebayes_1990_unique_snps, samtools_1990_unique_snps, gatk_1990_unique_snps, freebayes_2015_unique_snps, samtools_2015_unique_snps, gatk_2015_unique_snps)

norm_indels <- c(total_int_1990_indel, total_int_2015_indel, freebayes_1990_unique_indels, samtools_1990_unique_indels, gatk_1990_unique_indels, freebayes_2015_unique_indels, samtools_2015_unique_indels, gatk_2015_unique_indels)

#Calculate percents of snps and indels overlapping genes or repeats
percent_snps_gene <- as.data.frame(((colSums(snp_anno_gene_df[,4:ncol(snp_anno_gene_df)]))/norm_snps)*100)
percent_snps_gene <- cbind(rownames(percent_snps_gene), data.frame(percent_snps_gene, row.names=NULL))
colnames(percent_snps_gene) <- c("id", "Percent")
percent_snps_gene$Overlap <- "gene"

percent_indels_gene <- as.data.frame(((colSums(indel_anno_gene_df[,4:ncol(indel_anno_gene_df)]))/norm_indels)*100)
percent_indels_gene <- cbind(rownames(percent_indels_gene), data.frame(percent_indels_gene, row.names=NULL))
colnames(percent_indels_gene) <- c("id", "Percent")
percent_indels_gene$Overlap <- "gene"

percent_snps_repeat <- as.data.frame(((colSums(snp_anno_repeat_df[,4:ncol(snp_anno_repeat_df)]))/norm_snps)*100)
percent_snps_repeat <- cbind(rownames(percent_snps_repeat), data.frame(percent_snps_repeat, row.names=NULL))
colnames(percent_snps_repeat) <- c("id", "Percent")
percent_snps_repeat$Overlap <- "repeat"

percent_indels_repeat <- as.data.frame(((colSums(indel_anno_repeat_df[,4:ncol(indel_anno_repeat_df)]))/norm_indels)*100)
percent_indels_repeat <- cbind(rownames(percent_indels_repeat), data.frame(percent_indels_repeat, row.names=NULL))
colnames(percent_indels_repeat) <- c("id", "Percent")
percent_indels_repeat$Overlap <- "repeat"

#Do some reformatting of tables to get in better format for plotting
list <- list(percent_snps_gene, percent_snps_repeat, percent_indels_gene, percent_indels_repeat)
list <- lapply(list, mutate, Year=ifelse(grepl("1990",id), 1990, 2015))
list <- lapply(list, mutate, Variant=ifelse(grepl("SNP",id), "SNP", "INDEL"))
list <- lapply(list, function(x) {x$id = c("Shared", "Shared", "Unique FreeBayes", "Unique Samtools", "Unique GATK", "Unique FreeBayes", "Unique Samtools", "Unique GATK"); x})

#Combine all into dataframe based on shared column (id)
summary_table <- bind_rows(list)
summary_table$Year <- as.factor(summary_table$Year)

#Dot plots
#Set theme to be classic
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=12),
              axis.title = element_text(size=15),
              axis.line = element_line(size=0),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 13)
            )
)

plot <- ggplot(summary_table, aes(x = id, y = Percent)) +
  geom_point(aes(color = Overlap, shape = Year), alpha = 0.6, size = 3) +
  facet_wrap(~Variant) +
  xlab("") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = NA, color = "black", size = 1)) +
  scale_color_manual(values=c("darkgreen", "darkgoldenrod1"))
  

#Save the plot
pdf("variant_genomic_neighborhoods.pdf", width = 13)
plot
dev.off()
