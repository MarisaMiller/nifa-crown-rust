#UpSet diagram to compare overlap between variant callers
library(UpSetR)

#Set working directory to appropriate location
setwd("")

#Read in data
df <- read.delim("gatk_freebayes_samtools_compare_1990_2015.txt", header=TRUE)

#Pull out the relevant values for 1990 SNPs comparison
gatk_1990_unique_snps <- df[8,4] - df[11,4] - (df[12,4] - df[11,4]) - (df[13,4] -  df[11,4])
freebayes_1990_unique_snps <- df[9,4] - df[11,4] - (df[12,4] - df[11,4]) - (df[14,4] -  df[11,4])
samtools_1990_unique_snps <- df[10,4] - df[11,4] - (df[13,4] - df[11,4]) - (df[14,4] -  df[11,4])
intersect_G_FB_1990_snps <- df[12,4]
intersect_FB_ST_1990_snps <- df[14,4]
intersect_G_ST_1990_snps <- df[13,4]
intersect_all_1990_snps <- df[11,4]

snps_1990 <- c("GATK" = gatk_1990_unique_snps, "FreeBayes" = freebayes_1990_unique_snps, "Samtools" = samtools_1990_unique_snps, "GATK&FreeBayes" = intersect_G_FB_1990_snps, "FreeBayes&Samtools" = intersect_FB_ST_1990_snps, "GATK&Samtools" = intersect_G_ST_1990_snps, "GATK&FreeBayes&Samtools" = intersect_all_1990_snps)

upset(fromExpression(snps_1990), order.by = "freq", mainbar.y.label = "Variant Caller Intersections", sets.x.label = "Total Variants Per Caller", text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.1))
grid.edit('arrange', name='arrange2')
snps_1990_plot = grid.grab()

#Pull out the relevant values for 2015 SNPs comparison
gatk_2015_unique_snps <- df[1,4] - df[4,4] - (df[5,4] - df[4,4]) - (df[6,4] -  df[4,4])
freebayes_2015_unique_snps <- df[2,4] - df[4,4] - (df[5,4] - df[4,4]) - (df[7,4] -  df[4,4])
samtools_2015_unique_snps <- df[3,4] - df[4,4] - (df[6,4] - df[4,4]) - (df[7,4] -  df[4,4])
intersect_G_FB_2015_snps <- df[5,4]
intersect_FB_ST_2015_snps <- df[7,4]
intersect_G_ST_2015_snps <- df[6,4]
intersect_all_2015_snps <- df[4,4]

snps_2015 <- c("GATK" = gatk_2015_unique_snps, "FreeBayes" = freebayes_2015_unique_snps, "Samtools" = samtools_2015_unique_snps, "GATK&FreeBayes" = intersect_G_FB_2015_snps, "FreeBayes&Samtools" = intersect_FB_ST_2015_snps, "GATK&Samtools" = intersect_G_ST_2015_snps, "GATK&FreeBayes&Samtools" = intersect_all_2015_snps)

upset(fromExpression(snps_2015), order.by = "freq", mainbar.y.label = "Variant Caller Intersections", sets.x.label = "Total Variants Per Caller", text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.1))
grid.edit('arrange', name='arrange2')
snps_2015_plot = grid.grab()

#Pull out the relevant values for 1990 INDELs comparison
gatk_1990_unique_indels <- df[8,5] - df[11,5] - (df[12,5] - df[11,5]) - (df[13,5] -  df[11,5])
freebayes_1990_unique_indels <- df[9,5] - df[11,5] - (df[12,5] - df[11,5]) - (df[14,5] -  df[11,5])
samtools_1990_unique_indels <- df[10,5] - df[11,5] - (df[13,5] - df[11,5]) - (df[14,5] -  df[11,5])
intersect_G_FB_1990_indels <- df[12,5]
intersect_FB_ST_1990_indels <- df[14,5]
intersect_G_ST_1990_indels <- df[13,5]
intersect_all_1990_indels <- df[11,5]

indels_1990 <- c("GATK" = gatk_1990_unique_indels, "FreeBayes" = freebayes_1990_unique_indels, "Samtools" = samtools_1990_unique_indels, "GATK&FreeBayes" = intersect_G_FB_1990_indels, "FreeBayes&Samtools" = intersect_FB_ST_1990_indels, "GATK&Samtools" = intersect_G_ST_1990_indels, "GATK&FreeBayes&Samtools" = intersect_all_1990_indels)

upset(fromExpression(indels_1990), order.by = "freq", mainbar.y.label = "Variant Caller Intersections", sets.x.label = "Total Variants Per Caller", text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.1))
grid.edit('arrange', name='arrange2')
indels_1990_plot = grid.grab()

#Pull out the relevant values for 2015 indels comparison
gatk_2015_unique_indels <- df[1,5] - df[4,5] - (df[5,5] - df[4,5]) - (df[6,5] -  df[4,5])
freebayes_2015_unique_indels <- df[2,5] - df[4,5] - (df[5,5] - df[4,5]) - (df[7,5] -  df[4,5])
samtools_2015_unique_indels <- df[3,5] - df[4,5] - (df[6,5] - df[4,5]) - (df[7,5] -  df[4,5])
intersect_G_FB_2015_indels <- df[5,5]
intersect_FB_ST_2015_indels <- df[7,5]
intersect_G_ST_2015_indels <- df[6,5]
intersect_all_2015_indels <- df[4,5]

indels_2015 <- c("GATK" = gatk_2015_unique_indels, "FreeBayes" = freebayes_2015_unique_indels, "Samtools" = samtools_2015_unique_indels, "GATK&FreeBayes" = intersect_G_FB_2015_indels, "FreeBayes&Samtools" = intersect_FB_ST_2015_indels, "GATK&Samtools" = intersect_G_ST_2015_indels, "GATK&FreeBayes&Samtools" = intersect_all_2015_indels)

upset(fromExpression(indels_2015), order.by = "freq", mainbar.y.label = "Variant Caller Intersections", sets.x.label = "Total Variants Per Caller", text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.1))
grid.edit('arrange', name='arrange2')
indels_2015_plot = grid.grab()

#Save the plot
pdf("gatk_freebayes_samtools_comparison.pdf", height = 10, width = 15)
grid.arrange(arrangeGrob(snps_1990_plot, top = "SNPs 1990 Isolates"), arrangeGrob(snps_2015_plot, top = "SNPs 2015 Isolates"), arrangeGrob(indels_1990_plot, top = "INDELs 1990 Isolates"), arrangeGrob(indels_2015_plot, top = "INDELs 2015 Isolates"), ncol=2, nrow=2)
dev.off()

###Compare variants called with FreeBayes between 1990 and 2015, using more stringent filtering params for downstream analysis

#Read in data
df2 <- read.delim("1990_2015_freebayes_compare.txt", header=TRUE)

#Pull out the relevant values for snps
fb_2015_unique_snps <- df2[1,4] - df2[3,4]
fb_1990_unique_snps <- df2[2,4] - df2[3,4]
intersect_1990_2015_snps <- df2[3,4]

snps <- c("1990" = fb_1990_unique_snps, "2015" = fb_2015_unique_snps, "1990&2015" = intersect_1990_2015_snps)

snps_plot <- grid.grabExpr(upset(fromExpression(snps), order.by = "freq", mainbar.y.label = "Variant Intersections", sets.x.label = "Total Variants Per Year", text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.1)))

#Pull out the relevant values for indels
fb_2015_unique_indels <- df2[1,5] - df2[3,5]
fb_1990_unique_indels <- df2[2,5] - df2[3,5]
intersect_1990_2015_indels <- df2[3,5]

indels <- c("1990" = fb_1990_unique_indels, "2015" = fb_2015_unique_indels, "1990&2015" = intersect_1990_2015_indels)

indels_plot <- grid.grabExpr(upset(fromExpression(indels), order.by = "freq", mainbar.y.label = "Variant Intersections", sets.x.label = "Total Variants Per Year", text.scale = c(1.5, 1.5, 1, 1, 1.5, 1.1)))

#Save the plot
pdf("1990_2015_fb_comparison.pdf", height = 10, width = 6.5)
grid.arrange(arrangeGrob(snps_plot, top = "SNPs"), arrangeGrob(indels_plot, top = "INDELs"), nrow=2)
dev.off()
