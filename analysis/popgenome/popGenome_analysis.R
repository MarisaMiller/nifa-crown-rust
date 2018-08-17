#This script is for analysis with the PopGenome package
#Analysis performed with R 3.4.0 and PopGenome 2.2.4

#Change working directory to appropriate location
setwd("")

#Load required libraries
library(PopGenome)
library(ggplot2)
#library(ggridges) #Only needed if want to make ridgeline plots instead
library(gridExtra)
library(stringr)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)

#Read in vcf and gff files by providing folder name where the files are located.
#Make sure the file extensions are gff (not gff3) or else they won't be successfully loaded.
#First, open the source code for the readData function, and comment out the lines of code, and then hit save in the window that pops up.
#    #if (format == "VCF") {
#    SNP.DATA = TRUE
#    FAST = TRUE
#}
#if (SNP.DATA) {
#    big.data = TRUE
#}
#This ensures that the ff package is not employed, which seems to cause numerous issues with downstream analyses nad moving Rdata objects back and forth between machines
trace("readData", edit=TRUE)
vcf_1990 <- readData("1990/", format="VCF", gffpath="GFF_1990/", FAST = TRUE) #There are a few less GFF files for 1990 since there aren't VCF files for a couple of contigs
vcf_2015 <- readData("2015/", format="VCF", gffpath="GFF/", FAST = TRUE)
vcf_1990_2015 <- readData("1990_2015/", format="VCF", gffpath="GFF/", FAST = TRUE)

#Check some information about GENOME.class objects to make sure everything is set correctly
#get.sum.data(vcf_1990)
#get.sum.data(vcf_2015)

#Save and re-load objects for later analysis
#save(vcf_1990, file="vcf_1990.Rdata")
#save(vcf_2015, file="vcf_2015.Rdata")
#save(vcf_1990_2015, file="vcf_1990_2015.Rdata")
load("vcf_1990.Rdata")
load("vcf_2015.Rdata")
load("vcf_1990_2015.Rdata")

#Various functions for grabbing out values of analyses. I could probably make one function to do all of this, and just have more options, but this is good enough for now.
#Function to make dataframe for Nucleotide Diversity
#The values have to be normalized by the number of nucleotides in each window/region (could be whole contig level, gene length, or window size)
nucDiv.df <- function(x, year, gene_type=NULL) {
  if (is.null(gene_type))
    df <- data.frame("Nucleotide Diversity" = (x@nuc.diversity.within/x@n.sites), "Year" = year, "VCF" = x@region.names)
  else
    df <- data.frame("Nucleotide Diversity" = (x@nuc.diversity.within/x@n.sites), "Year" = year, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", x@region.names, perl=TRUE), "Gene.Type" = gene_type)
  df
}

#Function to make dataframe for Watterson's Theta
thetaW.df <- function(x, year, gene_type=NULL) {
  if (is.null(gene_type))
    df <- data.frame("Watterson's Theta" = (x@theta_Watterson/x@n.sites), "Year" = year, "VCF" = x@region.names)
  else
    df <- data.frame("Watterson's Theta" = (x@theta_Watterson/x@n.sites), "Year" = year, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", x@region.names, perl=TRUE), "Gene.Type" = gene_type)
  df
}

#Function to make dataframe for Tajima's D
tajD.df <- function(x, year, gene_type=NULL) {
  if (is.null(gene_type))
    df <- data.frame("Tajima's D" = x@Tajima.D, "Year" = year, "VCF" = x@region.names)
  else
    df <- data.frame("Tajima's D" = x@Tajima.D, "Year" = year, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", x@region.names, perl=TRUE), "Gene.Type" = gene_type)
  df
}

#Function to make dataframe for Fst
#nuc.F_ST.vs.all is the fixation index for each population against all individuals in other populations
fst.df <- function(x, year, gene_type=NULL) {
  if (is.null(gene_type))
    df <- data.frame("Fst" = x@nuc.F_ST.vs.all, "Year" = year, "VCF" = x@region.names)
  else
    df <- data.frame("Fst" = x@nuc.F_ST.vs.all, "Year" = year, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", x@region.names, perl=TRUE), "Gene.Type" = gene_type)
  df
}

#Plot settings
#Set theme to be classic for rest of plots
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=12),
              axis.title = element_text(size=15),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 13)
            )
)

######Genome-wide analysis for nucleotide diversity and Tajima's D######

#First analyze two years for diversity stats on genome-wide level
vcf_1990 <- diversity.stats(vcf_1990)
vcf_2015 <- diversity.stats(vcf_2015)

#Use nucDiv.df function to extract values
vcf_1990_nucdiv_norm <- nucDiv.df(vcf_1990, "1990")
vcf_2015_nucdiv_norm <- nucDiv.df(vcf_2015, "2015")

#Make a dataframe for plotting
nucdiv_df <- rbind(vcf_1990_nucdiv_norm, vcf_2015_nucdiv_norm)
colnames(nucdiv_df)[[1]] <- "Nucleotide Diversity" #Somehow I needed to rename this column again

#Get average across genome for 2 years
nuc_div_1990 <- mean(vcf_1990_nucdiv_norm$pop.1, na.rm = TRUE) #0.002550695
nuc_div_2015 <- mean(vcf_2015_nucdiv_norm$pop.1, na.rm = TRUE) #0.002154375

#Test if significant difference between two years using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pv_nuc_div <- wilcox.test(vcf_1990_nucdiv_norm$pop.1, vcf_2015_nucdiv_norm$pop.1) #p = 6.522e-08

#Analyze neutrality stats on genome-wide level
vcf_1990 <- neutrality.stats(vcf_1990)
vcf_2015 <- neutrality.stats(vcf_2015)

#Use thetaW.df function to extract values
vcf_1990_theta_norm <- thetaW.df(vcf_1990, "1990")
vcf_2015_theta_norm <- thetaW.df(vcf_2015, "2015")

vcf_1990_TajD <- tajD.df(vcf_1990, "1990")
vcf_2015_TajD <- tajD.df(vcf_2015, "2015")

#Make a dataframe for plotting
thetaW_df <- rbind(vcf_1990_theta_norm, vcf_2015_theta_norm)
colnames(thetaW_df)[[1]] <- "Watterson's Theta" #Somehow I needed to rename this column again

TajD_df <- rbind(vcf_1990_TajD, vcf_2015_TajD)
colnames(TajD_df)[[1]] <- "Tajima's D" #Somehow I needed to rename this column again

#Get average across genome for 2 years
thetaW_1990 <- mean(vcf_1990_theta_norm$pop.1, na.rm = TRUE) #0.001947765
thetaW_2015 <- mean(vcf_2015_theta_norm$pop.1, na.rm = TRUE) #0.00158565

TajD_1990 <- mean(vcf_1990_TajD$pop.1, na.rm = TRUE) #1.028332
TajD_2015 <- mean(vcf_2015_TajD$pop.1, na.rm = TRUE) #1.289528

#Test if significant difference between two years using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pv_thetaW <- wilcox.test(vcf_1990_theta_norm$pop.1, vcf_2015_theta_norm$pop.1) #p = 4.689e-13
pv_TajD <- wilcox.test(vcf_1990_TajD$pop.1, vcf_2015_TajD$pop.1) #p = 7.341e-08

#Test significant departure from neutrality by running Hudson’s MS program with the MS PopGenome module
#If no additional parameters are specified for the MS simulations, PopGenome will use the standard neutral model (SNM).
vcf_1990_coalesce <- MS(vcf_1990, thetaID = "Tajima", neutrality = TRUE, niter = 100)
vcf_2015_coalesce <- MS(vcf_2015, thetaID = "Tajima", neutrality = TRUE, niter = 100)

#Save and re-load objects for later analysis as the coalescence takes a while to run
#save(vcf_1990_coalesce, file="vcf_1990_coalesce.Rdata")
#save(vcf_2015_coalesce, file="vcf_2015_coalesce.Rdata")
load("vcf_1990_coalesce.Rdata")
load("vcf_2015_coalesce.Rdata")

#To examine if the simulated data is equal to the observed data (and hence no evidence for non-neutral evolution), examine the slot @prob.equal, and examine the average and variance at the end of the vector
#Our observed values are equal to the simulated values showing that there is no non-neutral evolution in the Pca populations
vcf_1990_coalesce@prob.less[[1]][,1] #Average p = 0.82728157, variance = 0.04491854, so obs !< sim
vcf_1990_coalesce@prob.equal[[1]][,1] #Average p = 0, variance = 0, so obs = sim
vcf_2015_coalesce@prob.less[[1]][,1] #Average p = 0.85654421, variance = 0.04387858, so obs !< sim
vcf_2015_coalesce@prob.equal[[1]][,1] #Average p =  2.919111e-05, variance = 5.121248e-07, so obs = sim

#Plots
viol_nucdiv_plot <- ggplot(nucdiv_df, aes(nucdiv_df$Year, nucdiv_df$`Nucleotide Diversity`)) +
  geom_violin(fill="lightblue", draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "Year", y="Nucleotide Diversity")

viol_theta_plot <- ggplot(thetaW_df, aes(thetaW_df$Year, thetaW_df$`Watterson's Theta`)) +
  geom_violin(fill="lightgreen", draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "Year", y="Watterson's Theta")

viol_TajD_plot <- ggplot(TajD_df, aes(TajD_df$Year,TajD_df$`Tajima's D`)) +
  geom_violin(fill="lightgoldenrod", draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "Year", y="Tajima's D")

#Alternative Plots using ggridges
# nucdiv_plot <- ggplot(nucdiv_df, aes(nucdiv_df$`Nucleotide Diversity`, nucdiv_df$Year, fill=factor(..quantile..))) +
#   geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, scale=0.8) +
#   scale_fill_brewer(palette = "Blues", guide=FALSE) + #Colors represent quartiles
#   scale_y_discrete(expand=c(0.001, 0)) + 
#   theme(strip.text.x = element_blank()) +
#   labs(x = "Nucleotide Diversity", y="Year")
# 
# theta_plot <- ggplot(thetaW_df, aes(thetaW_df$`Watterson's Theta`, thetaW_df$Year, fill=factor(..quantile..))) +
#   geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, scale=0.8) +
#   scale_fill_brewer(palette = "Greens", guide=FALSE) + #Colors represent quartiles
#   scale_y_discrete(expand=c(0.001, 0)) + 
#   theme(strip.text.x = element_blank()) +
#   labs(x = "Watterson's Theta", y="Year")
# 
# TajD_plot <- ggplot(TajD_df, aes(TajD_df$`Tajima's D`, TajD_df$Year, fill=factor(..quantile..))) +
#   geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, scale=0.8) +
#   scale_fill_brewer(palette = "Oranges", guide=FALSE) + #Colors represent quartiles
#   scale_y_discrete(expand=c(0.001, 0)) + 
#   theme(strip.text.x = element_blank()) +
#   labs(x = "Tajima's D", y="Year")

#Save the plots
pdf("viol_nuc_div_theta_TajD.pdf", height = 9, width = 5)
grid.arrange(arrangeGrob(viol_nucdiv_plot), arrangeGrob(viol_theta_plot), arrangeGrob(viol_TajD_plot), ncol=1, nrow=3)
dev.off()

# pdf("nuc_div_theta_TajD.pdf", height = 9, width = 5)
# grid.arrange(arrangeGrob(nucdiv_plot), arrangeGrob(theta_plot), arrangeGrob(TajD_plot), ncol=1, nrow=3)
# dev.off()




######Nucleotide diversity and Tajima's D for effectors vs all genes######

all_genes_1990 <- splitting.data(vcf_1990, subsites="gene", whole.data = FALSE)
all_genes_2015 <- splitting.data(vcf_2015, subsites="gene", whole.data = FALSE)

all_genes_1990 <- diversity.stats(all_genes_1990)
all_genes_2015 <- diversity.stats(all_genes_2015)

all_genes_1990 <- neutrality.stats(all_genes_1990)
all_genes_2015 <- neutrality.stats(all_genes_2015)

#Use function with gene_type argument
all_genes_1990_nucdiv_norm <- nucDiv.df(all_genes_1990, "1990", "All Genes")
all_genes_2015_nucdiv_norm <- nucDiv.df(all_genes_2015, "2015", "All Genes")

effector_gff_table   <- read.table("./gene_lists/effector_genes_sd80_p.gff3", sep="\t", colClasses=c("character","NULL","NULL","numeric","numeric",rep("NULL",3),"character"))

effector_parsing_list <- paste(effector_gff_table$V1, effector_gff_table$V4, effector_gff_table$V5, sep = ":")

effectors_1990_nucdiv_norm <- all_genes_1990_nucdiv_norm[all_genes_1990_nucdiv_norm$Contig.Start.Stop %in% effector_parsing_list, ]
effectors_2015_nucdiv_norm <- all_genes_2015_nucdiv_norm[all_genes_2015_nucdiv_norm$Contig.Start.Stop %in% effector_parsing_list, ]

effectors_1990_nucdiv_norm$Gene.Type <- "Effectors"
effectors_2015_nucdiv_norm$Gene.Type <- "Effectors"

#Some optional analysis with just haustorial effectors
# haus_effector_gff_table <- read.table("./gene_lists/haus_expr_secreted_effectors_sd80_p.gff3", sep="\t", colClasses=c("character","NULL","NULL","numeric","numeric",rep("NULL",3),"character"))
# 
# haus_effector_parsing_list <- paste(haus_effector_gff_table$V1, haus_effector_gff_table$V4, haus_effector_gff_table$V5, sep = ":")
# 
# haus_effectors_1990_nucdiv_norm <- all_genes_1990_nucdiv_norm[all_genes_1990_nucdiv_norm$Contig.Start.Stop %in% haus_effector_parsing_list, ]
# haus_effectors_2015_nucdiv_norm <- all_genes_2015_nucdiv_norm[all_genes_2015_nucdiv_norm$Contig.Start.Stop %in% haus_effector_parsing_list, ]
# 
# haus_effectors_1990_nucdiv_norm$Gene.Type <- "Haustorial Effectors"
# haus_effectors_2015_nucdiv_norm$Gene.Type <- "Haustorial Effectors"

#Make a dataframe for plotting
genes_effectors_nucdiv_df <- rbind(all_genes_1990_nucdiv_norm, effectors_1990_nucdiv_norm, all_genes_2015_nucdiv_norm, effectors_2015_nucdiv_norm)
# genes_effectors_nucdiv_df <- rbind(all_genes_1990_nucdiv_norm, effectors_1990_nucdiv_norm, haus_effectors_1990_nucdiv_norm, all_genes_2015_nucdiv_norm, effectors_2015_nucdiv_norm, haus_effectors_2015_nucdiv_norm)
colnames(genes_effectors_nucdiv_df)[[1]] <- "Nucleotide Diversity" #Somehow I needed to rename this column again

#Make a dataframe for effectors to export
effectors_nucdiv_df <- rbind(effectors_1990_nucdiv_norm, effectors_2015_nucdiv_norm)
colnames(effectors_nucdiv_df)[[1]] <- "Nucleotide Diversity"
effectors_nucdiv_df <- effectors_nucdiv_df %>% select(-Gene.Type)
#Format contig:coordinate format back to gene_ID and get format suitable for summary dataframe
effector_replace_list <- as.data.frame(paste(effector_gff_table$V1, effector_gff_table$V4, effector_gff_table$V5, sep = ":"))
colnames(effector_replace_list)[[1]] <- "Contig.Start.Stop"
effector_replace_list$gene_ID <- effector_gff_table$V9
effectors_nucdiv_df <- left_join(effector_replace_list, effectors_nucdiv_df, by = c("Contig.Start.Stop"))
wide_effectors_nucdiv_df <- spread(effectors_nucdiv_df, Year, "Nucleotide Diversity")
colnames(wide_effectors_nucdiv_df)[3:4] <- c("Nucleotide Diversity 1990", "Nucleotide Diversity 2015")

write.table(wide_effectors_nucdiv_df, file = "nucdiv_effectors.txt", sep = "\t", row.names = FALSE)

#Get average for effectors vs all genes for each year
all_gene_nuc_div_1990 <- mean(all_genes_1990_nucdiv_norm$pop.1, na.rm = TRUE) #0.002754378
effector_nuc_div_1990 <- mean(effectors_1990_nucdiv_norm$pop.1, na.rm = TRUE) #0.002915493
# haus_effector_nuc_div_1990 <- mean(haus_effectors_1990_nucdiv_norm$pop.1, na.rm = TRUE) #0.002411828

all_gene_nuc_div_2015 <- mean(all_genes_2015_nucdiv_norm$pop.1, na.rm = TRUE) #0.00236181
effector_nuc_div_2015 <- mean(effectors_2015_nucdiv_norm$pop.1, na.rm = TRUE) #0.002294652
# haus_effector_nuc_div_2015 <- mean(haus_effectors_2015_nucdiv_norm$pop.1, na.rm = TRUE) #0.001783428

#Test if significant difference between all genes and effectors using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pv_gene_1990_nuc_div <- wilcox.test(all_genes_1990_nucdiv_norm$pop.1, effectors_1990_nucdiv_norm$pop.1) #p = 0.0499
pv_gene_2015_nuc_div <- wilcox.test(all_genes_2015_nucdiv_norm$pop.1, effectors_2015_nucdiv_norm$pop.1) #p = 0.9763

# pv_gene_1990_nuc_div <- data.frame("All Genes vs Effectors" = wilcox.test(all_genes_1990_nucdiv_norm$pop.1, effectors_1990_nucdiv_norm$pop.1)$p.value, "All Genes vs Haustorial Effectors" = wilcox.test(all_genes_1990_nucdiv_norm$pop.1, haus_effectors_1990_nucdiv_norm$pop.1)$p.value, "Effectors vs Haustorial Effectors" = wilcox.test(effectors_1990_nucdiv_norm$pop.1, haus_effectors_1990_nucdiv_norm$pop.1)$p.value)
# All.Genes.vs.Effectors All.Genes.vs.Haustorial.Effectors Effectors.vs.Haustorial.Effectors
# 0.04990253                         0.9139987                         0.4537734 #sample size of haus effectors is too small to detect any significant different most likely

# pv_gene_2015_nuc_div <- data.frame("All Genes vs Effectors" = wilcox.test(all_genes_2015_nucdiv_norm$pop.1, effectors_2015_nucdiv_norm$pop.1)$p.value, "All Genes vs Haustorial Effectors" = wilcox.test(all_genes_2015_nucdiv_norm$pop.1, haus_effectors_2015_nucdiv_norm$pop.1)$p.value, "Effectors vs Haustorial Effectors" = wilcox.test(effectors_2015_nucdiv_norm$pop.1, haus_effectors_2015_nucdiv_norm$pop.1)$p.value)
# All.Genes.vs.Effectors All.Genes.vs.Haustorial.Effectors Effectors.vs.Haustorial.Effectors
# 0.9762626                         0.3748384                         0.4177304

#Neutrality stats
all_genes_1990_theta_norm <- thetaW.df(all_genes_1990, "1990", "All Genes")
all_genes_2015_theta_norm <- thetaW.df(all_genes_2015, "2015", "All Genes")

effectors_1990_theta_norm <- all_genes_1990_theta_norm[all_genes_1990_theta_norm$Contig.Start.Stop %in% effector_parsing_list, ]
effectors_2015_theta_norm <- all_genes_2015_theta_norm[all_genes_2015_theta_norm$Contig.Start.Stop %in% effector_parsing_list, ]

effectors_1990_theta_norm$Gene.Type <- "Effectors"
effectors_2015_theta_norm$Gene.Type <- "Effectors"

all_genes_1990_TajD <- tajD.df(all_genes_1990, "1990", "All Genes")
all_genes_2015_TajD <- tajD.df(all_genes_2015, "2015", "All Genes")

effectors_1990_TajD <- all_genes_1990_TajD[all_genes_1990_TajD$Contig.Start.Stop %in% effector_parsing_list, ]
effectors_2015_TajD <- all_genes_2015_TajD[all_genes_2015_TajD$Contig.Start.Stop %in% effector_parsing_list, ]

effectors_1990_TajD$Gene.Type <- "Effectors"
effectors_2015_TajD$Gene.Type <- "Effectors"

#Make dataframes for plotting
genes_effectors_thetaW_df <- rbind(all_genes_1990_theta_norm, effectors_1990_theta_norm, all_genes_2015_theta_norm, effectors_2015_theta_norm)
colnames(genes_effectors_thetaW_df)[[1]] <- "Watterson's Theta" #Somehow I needed to rename this column again

genes_effectors_TajD_df <- rbind(all_genes_1990_TajD, effectors_1990_TajD, all_genes_2015_TajD, effectors_2015_TajD)
colnames(genes_effectors_TajD_df)[[1]] <- "Tajima's D" #Somehow I needed to rename this column again

#Make dataframes for effectors to export
effectors_theta_df <- rbind(effectors_1990_theta_norm, effectors_2015_theta_norm)
colnames(effectors_theta_df)[[1]] <- "Watterson's Theta"
effectors_theta_df <- effectors_theta_df %>% select(-Gene.Type)
#Format contig:coordinate format back to gene_ID and get format suitable for summary dataframe
effectors_theta_df <- left_join(effector_replace_list, effectors_theta_df, by = c("Contig.Start.Stop"))
wide_effectors_theta_df <- spread(effectors_theta_df, Year, "Watterson's Theta")
colnames(wide_effectors_theta_df)[3:4] <- c("Watterson's Theta 1990", "Watterson's Theta 2015")

write.table(wide_effectors_theta_df, file = "theta_effectors.txt", sep = "\t", row.names = FALSE)

effectors_TajD_df <- rbind(effectors_1990_TajD, effectors_2015_TajD)
colnames(effectors_TajD_df)[[1]] <- "Tajima's D"
effectors_TajD_df <- effectors_TajD_df %>% select(-Gene.Type)
#Format contig:coordinate format back to gene_ID and get format suitable for summary dataframe
effectors_TajD_df <- left_join(effector_replace_list, effectors_TajD_df, by = c("Contig.Start.Stop"))
wide_effectors_TajD_df <- spread(effectors_TajD_df, Year, "Tajima's D")
colnames(wide_effectors_TajD_df)[3:4] <- c("Tajima's D 1990", "Tajima's D 2015")

write.table(wide_effectors_TajD_df, file = "TajD_effectors.txt", sep = "\t", row.names = FALSE)

#Get average for effectors vs all genes for each year
all_gene_thetaW_1990 <- mean(all_genes_1990_theta_norm$pop.1, na.rm = TRUE) #0.002570439
all_gene_thetaW_2015 <- mean(all_genes_2015_theta_norm$pop.1, na.rm = TRUE) #0.002171463

effector_thetaW_1990 <- mean(effectors_1990_theta_norm$pop.1, na.rm = TRUE) #0.002694588
effector_thetaW_2015 <- mean(effectors_2015_theta_norm$pop.1, na.rm = TRUE) #0.002228631

all_gene_TajD_1990 <- mean(all_genes_1990_TajD$pop.1, na.rm = TRUE) #0.7266541
all_gene_TajD_2015 <- mean(all_genes_2015_TajD$pop.1, na.rm = TRUE) #0.850953

effector_TajD_1990 <- mean(effectors_1990_TajD$pop.1, na.rm = TRUE) #0.6095681
effector_TajD_2015 <- mean(effectors_2015_TajD$pop.1, na.rm = TRUE) #0.6540236

#Test if significant difference between all genes and effectors using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pv_gene_1990_thetaW <- wilcox.test(all_genes_1990_theta_norm$pop.1, effectors_1990_theta_norm$pop.1) #p = 0.02266
pv_gene_2015_thetaW <- wilcox.test(all_genes_2015_theta_norm$pop.1, effectors_2015_theta_norm$pop.1) #p = 0.1238

pv_gene_1990_TajD <- wilcox.test(all_genes_1990_TajD$pop.1, effectors_1990_TajD$pop.1) #p = 0.03257
pv_gene_2015_TajD <- wilcox.test(all_genes_2015_TajD$pop.1, effectors_2015_TajD$pop.1) #p = 0.00274

#Make plots
viol_gene_nucdiv_plot <- ggplot(genes_effectors_nucdiv_df, aes(genes_effectors_nucdiv_df$Gene.Type, sqrt(genes_effectors_nucdiv_df$`Nucleotide Diversity` + 10^-10), fill = genes_effectors_nucdiv_df$Gene.Type)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("lightblue", alpha("lightblue", 0.5))) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.position = "none") +
  labs(x = "", y=expression(sqrt("Nucleotide Diversity + 10" ^ -10))) +
  facet_wrap(~Year)

viol_gene_theta_plot <- ggplot(genes_effectors_thetaW_df, aes(genes_effectors_thetaW_df$Gene.Type, sqrt(genes_effectors_thetaW_df$`Watterson's Theta` + 10^-10), fill = genes_effectors_thetaW_df$Gene.Type)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("lightgreen", alpha("lightgreen", 0.5))) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.position = "none") +
  labs(x = "",y=expression(sqrt("Watterson's Theta + 10" ^ -10))) +
  facet_wrap(~Year)

viol_gene_TajD_plot <- ggplot(genes_effectors_TajD_df, aes(genes_effectors_TajD_df$Gene.Type, genes_effectors_TajD_df$`Tajima's D`, fill = genes_effectors_TajD_df$Gene.Type)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("lightgoldenrod", alpha("lightgoldenrod", 0.5))) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.position = "none") +
  labs(x = "", y="Tajima's D") +
  facet_wrap(~Year)

# gene_theta_plot <- ggplot(genes_effectors_thetaW_df, aes(genes_effectors_thetaW_df$`Watterson's Theta`, genes_effectors_thetaW_df$Gene.Type, fill=factor(..quantile..))) +
#   geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, scale=4) + #no alpha allowed with gradient
#   scale_fill_brewer(palette = "Greens", guide=FALSE) + #Colors represent quartiles
#   scale_y_discrete(expand=c(0.001, 0)) +
#   #theme(strip.text.x = element_blank()) +
#   labs(x = "Watterson's Theta", y="") +
#   facet_wrap(~Year)
# 
# gene_TajD_plot <- ggplot(genes_effectors_TajD_df, aes(genes_effectors_TajD_df$`Tajima's D`, genes_effectors_TajD_df$Gene.Type, fill=factor(..quantile..))) +
#   geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, scale=4) + #no alpha allowed with gradient
#   scale_fill_brewer(palette = "Oranges", guide=FALSE) + #Colors represent quartiles
#   scale_y_discrete(expand=c(0.001, 0)) + 
#   #theme(strip.text.x = element_blank()) +
#   labs(x = "Tajima's D", y="") +
#   facet_wrap(~Year)
# 
# gene_nucdiv_plot <- ggplot(genes_effectors_nucdiv_df, aes(genes_effectors_nucdiv_df$`Nucleotide Diversity`, genes_effectors_nucdiv_df$Gene.Type, fill=factor(..quantile..))) +
#   geom_density_ridges_gradient(calc_ecdf = TRUE, quantiles = 4, scale=4) + #no alpha allowed with gradient
#   scale_fill_brewer(palette = "Blues", guide=FALSE) + #Colors represent quartiles
#   scale_y_discrete(expand=c(0.001, 0)) +
#   #theme(strip.text.x = element_blank()) +
#   labs(x = "Nucleotide Diversity", y="") +
#   facet_wrap(~Year)

#Save plots
pdf("viol_gene_nuc_div_theta_TajD.pdf", height = 9, width = 5)
grid.arrange(arrangeGrob(viol_gene_nucdiv_plot), arrangeGrob(viol_gene_theta_plot), arrangeGrob(viol_gene_TajD_plot), ncol=1, nrow=3)
dev.off()

# pdf("gene_nuc_div_theta_TajD.pdf", height = 9, width = 5)
# grid.arrange(arrangeGrob(gene_nucdiv_plot), arrangeGrob(gene_theta_plot), arrangeGrob(gene_TajD_plot), ncol=1, nrow=3)
# dev.off()






######Population based analysis######

#Set populations for Fst and additional analysis
no_bt_1990 <- c("90GA16-1","90TX52-1","90TX47-1","90KS101-1","90TX45-1","90LA38-1","90TX58-1","90AR100-1","90TX70-1")
bt_1990 <- c("90MN137-1","90MN148-1","90MN149-2","90MN152-1","90MN153-1","90WI131-1","90WI132-1","90SD164-1","90SD171-1","90SD172-1","90PA162-1")
mn_bt_1990 <- c("90MN1B-1","90MN2B-1","90MN3B-1","90MN5B-1","90MN7B-1","90MN8B-2","90MN9B-4","90MN13B-3","90MN14B-1","90MN17B-1")

no_bt_2015 <- c("15MS7-1","15TX3-1","15FL1-2","15Fl1-4")
bt_2015 <- c("15ND19-2","15ND19-5","15ND20-3","15ND20-4","15MN10-4","15MN10-5","15MN27-3","15MN23-1","15MN24-1","15MN25-3","15SD30-1","15SD30-3","15SD11-1","15SD11-2","15NE8-4","15NE8-5","15NE9-1","15NE9-3","15OH12-3")
mn_bt_2015 <- c("15MN13-4","15MN14-4","15MN15-3","15MN16-3","15MN17-5","15MN18-1","15MN18-3")

vcf_1990_pop <- set.populations(vcf_1990, list("NO BT"=no_bt_1990, "BT"=bt_1990, "MN BT"=mn_bt_1990), diploid=TRUE)
vcf_2015_pop <- set.populations(vcf_2015, list("NO BT"=no_bt_2015, "BT"=bt_2015, "MN BT"=mn_bt_2015), diploid=TRUE)

#Check some information about GENOME.class objects to make sure everything is set correctly
#get.individuals(vcf_1990_pop)
#get.individuals(vcf_2015_pop)
#vcf_1990_pop@populations
#vcf_2015_pop@populations

vcf_1990_pop <- diversity.stats(vcf_1990_pop)
vcf_2015_pop <- diversity.stats(vcf_2015_pop)

vcf_1990_pop <- neutrality.stats(vcf_1990_pop)
vcf_2015_pop <- neutrality.stats(vcf_2015_pop)

#Use nucDiv.df function
pop_1990_nucdiv_norm <- nucDiv.df(vcf_1990_pop, "1990")
colnames(pop_1990_nucdiv_norm)[1:3] <- c("NO BT", "BT", "MN BT")

pop_2015_nucdiv_norm <- nucDiv.df(vcf_2015_pop, "2015")
colnames(pop_2015_nucdiv_norm)[1:3] <- c("NO BT", "BT", "MN BT")

#Make a dataframe for plotting
pop_nucdiv_df <- rbind(pop_1990_nucdiv_norm, pop_2015_nucdiv_norm)
pop_nucdiv_df <- melt(pop_nucdiv_df, id.vars = c("Year", "VCF"), variable.name = "Population", value.name = "Nucleotide Diversity")
#I needed to re-order by Year because melt made the different years show up in alternating order for each population. This made the ggplot plots not display all the data.
pop_nucdiv_df <- with(pop_nucdiv_df, pop_nucdiv_df[order(Year),])

#Get average across genome for populations
pop_nuc_div_1990 <- data.frame("NO BT" = (mean(pop_1990_nucdiv_norm$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_1990_nucdiv_norm$BT, na.rm = TRUE)), "MN BT" = (mean(pop_1990_nucdiv_norm$`MN BT`, na.rm = TRUE)))
# NO.BT          BT       MN.BT
# 0.002581159 0.002486948 0.002517263
pop_nuc_div_2015 <- data.frame("NO BT" = (mean(pop_2015_nucdiv_norm$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_2015_nucdiv_norm$BT, na.rm = TRUE)), "MN BT" = (mean(pop_2015_nucdiv_norm$`MN BT`, na.rm = TRUE)))
# NO.BT        BT       MN.BT
# 0.002158965 0.002146401 0.002242132

#Test if significant difference between populations using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pop_1990_pv_nuc_div <- data.frame("NO BT vs BT" = wilcox.test(pop_1990_nucdiv_norm$`NO BT`, pop_1990_nucdiv_norm$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_1990_nucdiv_norm$`NO BT`, pop_1990_nucdiv_norm$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_1990_nucdiv_norm$`MN BT`, pop_1990_nucdiv_norm$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT MN.BT.vs.BT
# 0.206805      0.3858218   0.6985441

pop_2015_pv_nuc_div <- data.frame("NO BT vs BT" = wilcox.test(pop_2015_nucdiv_norm$`NO BT`, pop_2015_nucdiv_norm$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_2015_nucdiv_norm$`NO BT`, pop_2015_nucdiv_norm$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_2015_nucdiv_norm$`MN BT`, pop_2015_nucdiv_norm$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT MN.BT.vs.BT
# 0.8609917      0.1623113    0.198568

#Watterson's Theta
pop_1990_theta_norm <- thetaW.df(vcf_1990_pop, "1990")
colnames(pop_1990_theta_norm)[1:3] <- c("NO BT", "BT", "MN BT")

pop_2015_theta_norm <- thetaW.df(vcf_2015_pop, "2015")
colnames(pop_2015_theta_norm)[1:3] <- c("NO BT", "BT", "MN BT")

#Make a dataframe for plotting
pop_thetaW_df <- rbind(pop_1990_theta_norm, pop_2015_theta_norm)
pop_thetaW_df <- melt(pop_thetaW_df, id.vars = c("Year", "VCF"), variable.name = "Population", value.name = "Watterson's Theta")
pop_thetaW_df <- with(pop_thetaW_df, pop_thetaW_df[order(Year),])

#Get average across genome for populations
pop_thetaW_1990 <- data.frame("NO BT" = (mean(pop_1990_theta_norm$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_1990_theta_norm$BT, na.rm = TRUE)), "MN BT" = (mean(pop_1990_theta_norm$`MN BT`, na.rm = TRUE)))
# NO.BT          BT       MN.BT
# 0.002131589 0.002004334 0.002168642
pop_thetaW_2015 <- data.frame("NO BT" = (mean(pop_2015_theta_norm$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_2015_theta_norm$BT, na.rm = TRUE)), "MN BT" = (mean(pop_2015_theta_norm$`MN BT`, na.rm = TRUE)))
# NO.BT          BT     MN.BT
# 0.001900787 0.001619236 0.001944363

#Test if significant difference between populations using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pop_1990_pv_theta <- data.frame("NO BT vs BT" = wilcox.test(pop_1990_theta_norm$`NO BT`, pop_1990_theta_norm$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_1990_theta_norm$`NO BT`, pop_1990_theta_norm$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_1990_theta_norm$`MN BT`, pop_1990_theta_norm$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT MN.BT.vs.BT
# 0.04375712      0.6053646 0.009644018

pop_2015_pv_theta <- data.frame("NO BT vs BT" = wilcox.test(pop_2015_theta_norm$`NO BT`, pop_2015_theta_norm$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_2015_theta_norm$`NO BT`, pop_2015_theta_norm$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_2015_theta_norm$`MN BT`, pop_2015_theta_norm$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT  MN.BT.vs.BT
# 1.084197e-05      0.3221605 5.477268e-09

#Tajima's D
pop_1990_TajD <- tajD.df(vcf_1990_pop, "1990")
colnames(pop_1990_TajD)[1:3] <- c("NO BT", "BT", "MN BT")

pop_2015_TajD <- tajD.df(vcf_2015_pop, "2015")
colnames(pop_2015_TajD)[1:3] <- c("NO BT", "BT", "MN BT")

#Make a dataframe for plotting
pop_TajD_df <- rbind(pop_1990_TajD, pop_2015_TajD)
pop_TajD_df <- melt(pop_TajD_df, id.vars = c("Year", "VCF"), variable.name = "Population", value.name = "Tajima's D")
pop_TajD_df <- with(pop_TajD_df, pop_TajD_df[order(Year),])

#Get average across genome for populations
pop_TajD_1990 <- data.frame("NO BT" = (mean(pop_1990_TajD$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_1990_TajD$BT, na.rm = TRUE)), "MN BT" = (mean(pop_1990_TajD$`MN BT`, na.rm = TRUE)))
# NO.BT        BT    MN.BT
# 0.8192679 0.9044047 0.636858

pop_TajD_2015 <- data.frame("NO BT" = (mean(pop_2015_TajD$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_2015_TajD$BT, na.rm = TRUE)), "MN BT" = (mean(pop_2015_TajD$`MN BT`, na.rm = TRUE)))
# NO.BT       BT     MN.BT
# 0.7818182 1.228937 0.7114847

#Test if significant difference between populations using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pop_1990_pv_TajD <- data.frame("NO BT vs BT" = wilcox.test(pop_1990_TajD$`NO BT`, pop_1990_TajD$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_1990_TajD$`NO BT`, pop_1990_TajD$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_1990_TajD$`MN BT`, pop_1990_TajD$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT  MN.BT.vs.BT
# 0.06848163   5.745884e-09 4.583158e-14

pop_2015_pv_TajD <- data.frame("NO BT vs BT" = wilcox.test(pop_2015_TajD$`NO BT`, pop_2015_TajD$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_2015_TajD$`NO BT`, pop_2015_TajD$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_2015_TajD$`MN BT`, pop_2015_TajD$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT  MN.BT.vs.BT
# 5.310849e-28     0.08822465 8.397051e-32

#Make plots of above
pop_nucdiv_plot <- ggplot(pop_nucdiv_df, aes(pop_nucdiv_df$Population, pop_nucdiv_df$`Nucleotide Diversity`, fill = pop_nucdiv_df$Population)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "", y="Nucleotide Diversity") +
  facet_wrap(~Year)

pop_theta_plot <- ggplot(pop_thetaW_df, aes(pop_thetaW_df$Population, pop_thetaW_df$`Watterson's Theta`, fill = pop_thetaW_df$Population)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "", y="Watterson's Theta") +
  facet_wrap(~Year)

pop_TajD_plot <- ggplot(pop_TajD_df, aes(pop_TajD_df$Population, pop_TajD_df$`Tajima's D`, fill = pop_TajD_df$Population)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "", y="Tajima's D") +
  facet_wrap(~Year)

#Save plots
pdf("pop_nuc_div_theta_TajD.pdf", height = 9, width = 5)
grid.arrange(arrangeGrob(pop_nucdiv_plot), arrangeGrob(pop_theta_plot), arrangeGrob(pop_TajD_plot), ncol=1, nrow=3)
dev.off()






######Fst######

vcf_1990_pop <- F_ST.stats(vcf_1990_pop, mode="nucleotide")
vcf_2015_pop <- F_ST.stats(vcf_2015_pop, mode="nucleotide")

#Here is the set up for Fst analysis using the data combined and using the year information as the population
year_1990 <- c("90GA16-1","90TX52-1","90TX47-1","90KS101-1","90TX45-1","90LA38-1","90TX58-1","90AR100-1","90TX70-1", "90MN137-1","90MN148-1","90MN149-2","90MN152-1","90MN153-1","90WI131-1","90WI132-1","90SD164-1","90SD171-1","90SD172-1","90PA162-1", "90MN1B-1","90MN2B-1","90MN3B-1","90MN5B-1","90MN7B-1","90MN8B-2","90MN9B-4","90MN13B-3","90MN14B-1","90MN17B-1")

year_2015 <- c("15MS7-1","15TX3-1","15FL1-2","15Fl1-4", "15ND19-2","15ND19-5","15ND20-3","15ND20-4","15MN10-4","15MN10-5","15MN27-3","15MN23-1","15MN24-1","15MN25-3","15SD30-1","15SD30-3","15SD11-1","15SD11-2","15NE8-4","15NE8-5","15NE9-1","15NE9-3","15OH12-3", "15MN13-4","15MN14-4","15MN15-3","15MN16-3","15MN17-5","15MN18-1","15MN18-3")

vcf_1990_2015_pop <- set.populations(vcf_1990_2015, list("1990"=year_1990, "2015"=year_2015), diploid=TRUE)
vcf_1990_2015_pop <- F_ST.stats(vcf_1990_2015_pop, mode="nucleotide")

#Use fst.df function
pop_1990_Fst <- fst.df(vcf_1990_pop, "1990")
colnames(pop_1990_Fst)[1:3] <- c("NO BT", "BT", "MN BT")

pop_2015_Fst <- fst.df(vcf_2015_pop, "2015")
colnames(pop_2015_Fst)[1:3] <- c("NO BT", "BT", "MN BT")

pop_1990_2015_Fst <- fst.df(vcf_1990_2015_pop, "1990_2015")
colnames(pop_1990_2015_Fst)[1:2] <- c("1990", "2015")

#Make a dataframe for plotting
pop_Fst_df <- rbind(pop_1990_Fst, pop_2015_Fst)
pop_Fst_df <- melt(pop_Fst_df, id.vars = c("Year", "VCF"), variable.name = "Population", value.name = "Fst")
pop_Fst_df <- with(pop_Fst_df, pop_Fst_df[order(Year),])

pop_combined_Fst_df <- pop_1990_2015_Fst %>% select(`1990`) %>% rename("Fst" = `1990`)

#Get average across genome for populations
pop_Fst_1990 <- data.frame("NO BT" = (mean(pop_1990_Fst$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_1990_Fst$BT, na.rm = TRUE)), "MN BT" = (mean(pop_1990_Fst$`MN BT`, na.rm = TRUE)))
# NO.BT         BT      MN.BT
# 0.02307465 0.01781665 0.00976822

pop_Fst_2015 <- data.frame("NO BT" = (mean(pop_2015_Fst$`NO BT`, na.rm = TRUE)), "BT" = (mean(pop_2015_Fst$BT, na.rm = TRUE)), "MN BT" = (mean(pop_2015_Fst$`MN BT`, na.rm = TRUE)))
# NO.BT          BT       MN.BT
# -0.01883792 -0.01517778 -0.01029115

pop_Fst_1990_2015 <- data.frame("Fst" = (mean(pop_1990_2015_Fst$`1990`, na.rm = TRUE))) #Only have two populations so each column is identical, just use one column for mean
# Fst
# 0.03933639

#Test if significant difference between populations using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pop_1990_pv_Fst <- data.frame("NO BT vs BT" = wilcox.test(pop_1990_Fst$`NO BT`, pop_1990_Fst$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_1990_Fst$`NO BT`, pop_1990_Fst$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_1990_Fst$`MN BT`, pop_1990_Fst$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT  MN.BT.vs.BT
# 0.001888792   2.126675e-12 0.0001376838

pop_2015_pv_Fst <- data.frame("NO BT vs BT" = wilcox.test(pop_2015_Fst$`NO BT`, pop_2015_Fst$BT)$p.value, "NO BT vs MN BT" = wilcox.test(pop_2015_Fst$`NO BT`, pop_2015_Fst$`MN BT`)$p.value, "MN BT vs BT" = wilcox.test(pop_2015_Fst$`MN BT`, pop_2015_Fst$BT)$p.value)
# NO.BT.vs.BT NO.BT.vs.MN.BT MN.BT.vs.BT
# 71.032128e-09   1.077967e-09   0.5077578

pop_year_compare_pv_Fst <- data.frame("NO BT" = wilcox.test(pop_1990_Fst$`NO BT`, pop_2015_Fst$`NO BT`)$p.value, "MN BT" = wilcox.test(pop_1990_Fst$`MN BT`, pop_2015_Fst$`MN BT`)$p.value, "BT" = wilcox.test(pop_1990_Fst$BT, pop_2015_Fst$BT)$p.value)
# NO.BT        MN.BT           BT
# 5.72018e-83 1.808471e-46 3.391648e-80

#Make plots of above
viol_pop_Fst_plot <- ggplot(pop_Fst_df, aes(pop_Fst_df$Population, pop_Fst_df$Fst, fill = pop_Fst_df$Population)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "", y=expression("F"[st])) +
  facet_wrap(~Year)

viol_pop_combined_Fst_plot <- ggplot(pop_combined_Fst_df, aes("", Fst)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "", y=expression(italic("F")[ST]))

#Let's look at pairwise Fst also to make a table
pairwise_Fst_1990 <- colMeans(t(vcf_1990_pop@nuc.F_ST.pairwise), na.rm = TRUE)
#"NO BT vs BT" "NO BT vs MN BT" "BT vs MN BT"
#0.030682091 0.014772420 0.003520168

pairwise_Fst_2015 <- colMeans(t(vcf_2015_pop@nuc.F_ST.pairwise), na.rm = TRUE)
#"NO BT vs BT" "NO BT vs MN BT" "BT vs MN BT"
#-0.024752221 -0.015041716 -0.006338667

pairwise_Fst_1990_2015 <- colMeans(t(vcf_1990_2015_pop@nuc.F_ST.pairwise), na.rm = TRUE) #since only two populations this is identical to above

#Save plot
pdf("pop_Fst.pdf", height = 5, width = 5)
viol_pop_Fst_plot
dev.off()

pdf("pop_combined_Fst.pdf", height = 3, width = 3)
viol_pop_combined_Fst_plot
dev.off()

#Split populations into genes, and look at nucleotide diversity in different populations for different classes
pop_all_genes_1990 <- splitting.data(vcf_1990_pop, subsites="gene", whole.data = FALSE)
pop_all_genes_2015 <- splitting.data(vcf_2015_pop, subsites="gene", whole.data = FALSE)
pop_all_genes_1990_2015 <- splitting.data(vcf_1990_2015_pop, subsites="gene", whole.data = FALSE)

pop_all_genes_1990 <- diversity.stats(pop_all_genes_1990)
pop_all_genes_2015 <- diversity.stats(pop_all_genes_2015)

pop_all_genes_1990_nucdiv_norm <- nucDiv.df(pop_all_genes_1990, "1990", "All Genes")
colnames(pop_all_genes_1990_nucdiv_norm)[1:3] <- c("NO BT", "BT", "MN BT")

pop_all_genes_2015_nucdiv_norm <- nucDiv.df(pop_all_genes_2015, "2015", "All Genes")
colnames(pop_all_genes_2015_nucdiv_norm)[1:3] <- c("NO BT", "BT", "MN BT")

pop_effectors_1990_nucdiv_norm <- pop_all_genes_1990_nucdiv_norm[pop_all_genes_1990_nucdiv_norm$Contig.Start.Stop %in% effector_parsing_list, ]
pop_effectors_2015_nucdiv_norm <- pop_all_genes_2015_nucdiv_norm[pop_all_genes_2015_nucdiv_norm$Contig.Start.Stop %in% effector_parsing_list, ]

pop_effectors_1990_nucdiv_norm$Gene.Type <- "Effectors"
pop_effectors_2015_nucdiv_norm$Gene.Type <- "Effectors"

#Make a dataframe for plotting
pop_genes_effectors_nucdiv_df <- rbind(pop_all_genes_1990_nucdiv_norm, pop_effectors_1990_nucdiv_norm, pop_all_genes_2015_nucdiv_norm, pop_effectors_2015_nucdiv_norm)

pop_genes_effectors_nucdiv_df <- melt(pop_genes_effectors_nucdiv_df, id.vars = c("Year", "Contig.Start.Stop", "Gene.Type"), variable.name = "Population", value.name = "Nucleotide Diversity")

pop_genes_effectors_nucdiv_df <- with(pop_genes_effectors_nucdiv_df, pop_genes_effectors_nucdiv_df[order(Year),])

#Get average across genome for populations
pop_genes_nuc_div_1990 <- data.frame("All Genes NO BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], na.rm = TRUE)), "All Genes BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], na.rm = TRUE)), "All Genes MN BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], na.rm = TRUE)))
# All.Genes.NO.BT All.Genes.BT All.Genes.MN.BT
# 0.002801132  0.002675879     0.002724328

pop_effectors_nuc_div_1990 <- data.frame("Effectors NO BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], na.rm = TRUE)), "Effectors BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], na.rm = TRUE)), "Effectors MN BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], na.rm = TRUE)))
# Effectors.NO.BT Effectors.BT Effectors.MN.BT
# 0.002956242  0.002783779      0.00290573

pop_genes_nuc_div_2015 <- data.frame("All Genes NO BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], na.rm = TRUE)), "All Genes BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], na.rm = TRUE)), "All Genes MN BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], na.rm = TRUE)))
# All.Genes.NO.BT All.Genes.BT All.Genes.MN.BT
# 0.00235414  0.002351131     0.002467288

pop_effectors_nuc_div_2015 <- data.frame("Effectors NO BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], na.rm = TRUE)), "Effectors BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], na.rm = TRUE)), "Effectors MN BT" = (mean(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], na.rm = TRUE)))
# Effectors.NO.BT Effectors.BT Effectors.MN.BT
# 0.002227317  0.002288395     0.002391199

#Test if significant difference between genes in different populations using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pop_all_genes_1990_pv_nuc_div <- data.frame("NO BT" = wilcox.test(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"])$p.value,
                                            "BT" = wilcox.test(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"])$p.value,
                                            "MN BT" = wilcox.test(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"], pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 1990, "Nucleotide Diversity"])$p.value)
# NO.BT        BT      MN.BT
# 0.1092406 0.1390622 0.02918518

pop_all_genes_2015_pv_nuc_div <- data.frame("NO BT" = wilcox.test(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "NO BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"])$p.value,
                                            "BT" = wilcox.test(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"])$p.value,
                                            "MN BT" = wilcox.test(pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "All Genes" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"], pop_genes_effectors_nucdiv_df[pop_genes_effectors_nucdiv_df$Population == "MN BT" & pop_genes_effectors_nucdiv_df$Gene.Type == "Effectors" & pop_genes_effectors_nucdiv_df$Year == 2015, "Nucleotide Diversity"])$p.value)
# NO.BT        BT     MN.BT
# 0.121567 0.9520147 0.8865135

#Make plot of above, but this plot isn't particularly useful, just a table is adequate for these results
viol_pop_gene_nucdiv_plot <- ggplot(pop_genes_effectors_nucdiv_df, aes(pop_genes_effectors_nucdiv_df$Population, sqrt(pop_genes_effectors_nucdiv_df$`Nucleotide Diversity` + 10^-10), fill = pop_genes_effectors_nucdiv_df$Gene.Type)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = c("lightblue", alpha("lightblue", 0.5))) +
  stat_summary(aes(color=pop_genes_effectors_nucdiv_df$Gene.Type),fun.y = mean, geom="point", size=1, position=position_dodge(.9), color="black") +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  labs(x = "", y=expression(sqrt("Nucleotide Diversity + 10" ^ -10))) +
  facet_wrap(~Year)

#Save plot
pdf("viol_pop_gene_nucdiv.pdf", height = 5, width = 5)
viol_pop_gene_nucdiv_plot
dev.off()


#Look at Fst values on a per-gene basis in the populations
#Do effectors have highest Fst values?
pop_all_genes_1990 <- F_ST.stats(pop_all_genes_1990, mode="nucleotide")
pop_all_genes_2015 <- F_ST.stats(pop_all_genes_2015, mode="nucleotide")
pop_all_genes_1990_2015 <- F_ST.stats(pop_all_genes_1990_2015, mode="nucleotide")

pop_all_genes_1990_Fst <- fst.df(pop_all_genes_1990, "1990", "All Genes")
colnames(pop_all_genes_1990_Fst)[1:3] <- c("NO BT", "BT", "MN BT")

pop_all_genes_2015_Fst <- fst.df(pop_all_genes_2015, "2015", "All Genes")
colnames(pop_all_genes_2015_Fst)[1:3] <- c("NO BT", "BT", "MN BT")

pop_all_genes_1990_2015_Fst <- fst.df(pop_all_genes_1990_2015, "1990_2015", "All Genes")
colnames(pop_all_genes_1990_2015_Fst)[1:2] <- c("Fst", "2015")

pop_effectors_1990_Fst <- pop_all_genes_1990_Fst[pop_all_genes_1990_Fst$Contig.Start.Stop %in% effector_parsing_list, ]
pop_effectors_2015_Fst <- pop_all_genes_2015_Fst[pop_all_genes_2015_Fst$Contig.Start.Stop %in% effector_parsing_list, ]
pop_effectors_1990_2015_Fst <- pop_all_genes_1990_2015_Fst[pop_all_genes_1990_2015_Fst$Contig.Start.Stop %in% effector_parsing_list, ]

pop_effectors_1990_Fst$Gene.Type <- "Effectors"
pop_effectors_2015_Fst$Gene.Type <- "Effectors"
pop_effectors_1990_2015_Fst$Gene.Type <- "Effectors"

#Make a dataframe for plotting
pop_genes_effectors_Fst_df <- rbind(pop_all_genes_1990_Fst, pop_effectors_1990_Fst, pop_all_genes_2015_Fst, pop_effectors_2015_Fst)
pop_genes_effectors_Fst_df <- melt(pop_genes_effectors_Fst_df, id.vars = c("Year", "Contig.Start.Stop", "Gene.Type"), variable.name = "Population", value.name = "Fst")
pop_genes_effectors_Fst_df <- with(pop_genes_effectors_Fst_df, pop_genes_effectors_Fst_df[order(Year),])

#Make dataframes for effectors to export
effectors_Fst_df <- rbind(pop_effectors_1990_Fst, pop_effectors_2015_Fst)
colnames(effectors_Fst_df)[1:3] <- c("NO BT Fst", "BT Fst", "MN BT Fst")
effectors_Fst_df <- effectors_Fst_df %>% select(-Gene.Type)
#Format contig:coordinate format back to gene_ID and get format suitable for summary dataframe
wide_effectors_Fst_df <- dcast(setDT(effectors_Fst_df), Contig.Start.Stop ~ Year, value.var = c("NO BT Fst", "BT Fst", "MN BT Fst")) 
wide_effectors_Fst_df <- left_join(effector_replace_list, wide_effectors_Fst_df, by = c("Contig.Start.Stop"))

write.table(wide_effectors_Fst_df, file = "Fst_effectors.txt", sep = "\t", row.names = FALSE)

effectors_combined_Fst_df <- pop_effectors_1990_2015_Fst %>% select(-c(Gene.Type, `2015`))
#Format contig:coordinate format back to gene_ID and get format suitable for summary dataframe
effectors_combined_Fst_df <- left_join(effector_replace_list, effectors_combined_Fst_df, by = c("Contig.Start.Stop"))

write.table(effectors_combined_Fst_df, file = "Fst_combined_effectors.txt", sep = "\t", row.names = FALSE)

#Get average across genome for populations
pop_genes_Fst_1990 <- data.frame("All Genes NO BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], na.rm = TRUE)), "All Genes BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], na.rm = TRUE)), "All Genes MN BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], na.rm = TRUE)))
# All.Genes.NO.BT All.Genes.BT All.Genes.MN.BT
# 0.01884833   0.01438181     0.009584078

pop_effectors_Fst_1990 <- data.frame("Effectors NO BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], na.rm = TRUE)), "Effectors BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], na.rm = TRUE)), "Effectors MN BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], na.rm = TRUE)))
# Effectors.NO.BT Effectors.BT Effectors.MN.BT
# 0.02016218    0.0144455     0.009485931

pop_genes_Fst_2015 <- data.frame("All Genes NO BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], na.rm = TRUE)), "All Genes BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], na.rm = TRUE)), "All Genes MN BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], na.rm = TRUE)))
# All.Genes.NO.BT All.Genes.BT All.Genes.MN.BT
# -0.006253911 -0.006769668    -0.003722236

pop_effectors_Fst_2015 <- data.frame("Effectors NO BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], na.rm = TRUE)), "Effectors BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], na.rm = TRUE)), "Effectors MN BT" = (mean(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], na.rm = TRUE)))
# Effectors.NO.BT Effectors.BT Effectors.MN.BT
# -1.412354e-05 -0.004549002     0.001928567

#Test if significant difference between genes in different populations using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
pop_all_genes_1990_pv_Fst <- data.frame("NO BT" = wilcox.test(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"])$p.value,
                                            "BT" = wilcox.test(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"])$p.value,
                                            "MN BT" = wilcox.test(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"], pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 1990, "Fst"])$p.value)
# NO.BT        BT     MN.BT
# 0.7924743 0.5918064 0.4401975

pop_all_genes_2015_pv_Fst <- data.frame("NO BT" = wilcox.test(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "NO BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"])$p.value,
                                        "BT" = wilcox.test(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"])$p.value,
                                        "MN BT" = wilcox.test(pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "All Genes" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"], pop_genes_effectors_Fst_df[pop_genes_effectors_Fst_df$Population == "MN BT" & pop_genes_effectors_Fst_df$Gene.Type == "Effectors" & pop_genes_effectors_Fst_df$Year == 2015, "Fst"])$p.value)
# NO.BT        BT      MN.BT
# 0.01486869 0.1449332 0.01412985

#Make plot of above, although this plot is not very useful for the data above, and a table is the best choice.
viol_pop_gene_Fst_plot <- ggplot(pop_genes_effectors_Fst_df, aes(pop_genes_effectors_Fst_df$Population, pop_genes_effectors_Fst_df$Fst, fill = pop_genes_effectors_Fst_df$Gene.Type)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  #stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  stat_summary(aes(color=pop_genes_effectors_Fst_df$Gene.Type),fun.y = mean, geom="point", size=1, position=position_dodge(.9), color="black") +
  scale_fill_manual(values = c("white", "grey")) +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  #geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.5) +
  labs(x = "", y=expression("F"[st])) +
  facet_wrap(~Year)

#Save plot
pdf("viol_pop_gene_Fst.pdf", height = 7, width = 7)
viol_pop_gene_Fst_plot
dev.off()



#######Selective sweep tests##########
## CLR
vcf_1990 <- detail.stats(vcf_1990)
freq_1990 <- vcf_1990@region.stats@minor.allele.freqs[[1]]
freq.table_1990 <- list()
freq.table_1990[[1]] <- table(freq_1990)

vcf_2015 <- detail.stats(vcf_2015)
freq_2015 <- vcf_2015@region.stats@minor.allele.freqs[[1]]
freq.table_2015 <- list()
freq.table_2015[[1]] <- table(freq_2015)

## calculate CLR
all_genes_1990_sweep <- sweeps.stats(all_genes_1990, freq.table=freq.table_1990)
all_genes_2015_sweep <- sweeps.stats(all_genes_2015, freq.table=freq.table_2015)

#Make data frame
clr_1990 <- data.frame("CLR" = all_genes_1990_sweep@CLR, "Year" = 1990, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", all_genes_1990_sweep@region.names, perl=TRUE), "Gene.Type" = "All Genes")

clr_2015 <- data.frame("CLR" = all_genes_2015_sweep@CLR, "Year" = 2015, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", all_genes_2015_sweep@region.names, perl=TRUE), "Gene.Type" = "All Genes")

effectors_clr_1990 <- clr_1990[clr_1990$Contig.Start.Stop %in% effector_parsing_list, ]
effectors_clr_2015 <- clr_2015[clr_2015$Contig.Start.Stop %in% effector_parsing_list, ]

effectors_clr_1990$Gene.Type <- "Effectors"
effectors_clr_2015$Gene.Type <- "Effectors"

#Make a dataframe for effectors to export
effectors_clr_df <- rbind(effectors_clr_1990, effectors_clr_2015)
colnames(effectors_clr_df)[[1]] <- "CLR"
effectors_clr_df <- effectors_clr_df %>% select(-Gene.Type)
#Format contig:coordinate format back to gene_ID and get format suitable for summary dataframe
effector_replace_list <- as.data.frame(paste(effector_gff_table$V1, effector_gff_table$V4, effector_gff_table$V5, sep = ":"))
colnames(effector_replace_list)[[1]] <- "Contig.Start.Stop"
effector_replace_list$gene_ID <- effector_gff_table$V9
effectors_clr_df <- left_join(effector_replace_list, effectors_clr_df, by = c("Contig.Start.Stop"))
wide_effectors_clr_df <- spread(effectors_clr_df, Year, "CLR")
colnames(wide_effectors_clr_df)[3:4] <- c("CLR 1990", "CLR 2015")

write.table(wide_effectors_clr_df, file = "clr_effectors.txt", sep = "\t", row.names = FALSE)















#The sections below are for analysis ideas that I experimented with or that failed

#######Sliding Windows##########
#Generate sliding windows
#The options below will scan the data with 10,000 bp consecutive nucleotide windows.
#vcf_1990_slide <- sliding.window.transform(vcf_1990, 10000, 10000, type=2,  whole.data = FALSE) #FALSE can only be used if big.data = FALSE

#vcf_1990_slide <- diversity.stats(vcf_1990_slide)
#vcf_1990_slide@nuc.diversity.within #divide by window length to get nuc div per window

#Should I divide genome into sliding windows, and look at windows with effectors vs without effectors? Or, extract highest "outlier" windows (based on some arbitrary cutoff), and then see how many effectors are there?
#If interesting signal, I should plot a window across a candidate and see what it looks like. Perhaps an avr homolog? Not AvrSr50 since it has no homologs, but what about 35 or 27??

#######Synonymous and non-synonymous variants##########
#After further consideration, I'm not sure if using the variant annotation within PopGenome is the best approach. I decided to use Annovar for variant classification, and then used variants annotated with that method for DAPC analysis.
#Also, since I can't perform the MK test because it requires an outgroup, I don't need to classify syn/nonsyn variants within PopGenome.
#Below are some examples for performing syn/nonsyn classification with PopGenome that I tested prior to determining I will not be using these analyses.

#Set syn and non-syn variants for downstream analysis (if doing the MK test, etc) and for simple classification
#This has to be done on MSI as on my local computer the calculation stalls even if I use the ff package. However, if I don't use the ff package, syn and non-syn variants cannot be calculated because certain slots aren't filled unless big.data = TRUE.

#Because the vcf files for 1990 do not include variants for 2 of the contigs, a folder was created for just these FASTAs to prevent a "subscript out of bounds error" for 1990 isolates
#contig_files_1990 <- dput(as.character(list.files("../FASTA_1990", full.names = TRUE)))
#contig_files <- dput(as.character(list.files("../FASTA", full.names = TRUE)))

#vcf_1990_syn <- set.synnonsyn(vcf_1990, ref.chr = contig_files_1990)
#vcf_2015_syn <- set.synnonsyn(vcf_2015, ref.chr = contig_files)

#save(vcf_1990_syn, file="vcf_1990_syn.Rdata")
#save(vcf_2015_syn, file="vcf_2015_syn.Rdata")

#all_genes_1990_syn <- set.synnonsyn(all_genes_1990, ref.chr = contig_files_1990)
#all_genes_2015_syn <- set.synnonsyn(all_genes_2015, ref.chr = contig_files)

#save(all_genes_1990_syn, file="all_genes_1990_syn.Rdata")
#save(all_genes_2015_syn, file="all_genes_2015_syn.Rdata")

#Then just load the Rdata saved on MSI to summarize variants
# load("vcf_1990_syn.Rdata")
# load("vcf_2015_syn.Rdata")

# number of synonymous changes in whole genome as classified by PopGenome
# syn_list_1990 <- lapply(vcf_1990_syn@region.data@synonymous, function(x) sum(x==1, na.rm=TRUE))
# syn_total_1990 <- Reduce("+", syn_list_1990) #58,020 variants
# syn_list_2015 <- lapply(vcf_2015_syn@region.data@synonymous, function(x) sum(x==1, na.rm=TRUE))
# syn_total_2015 <- Reduce("+", syn_list_2015) #49,215 variants

# number of non-synonymous changes in whole genome
# nonsyn_list_1990 <- lapply(vcf_1990@region.data@synonymous, function(x) sum(x==0, na.rm=TRUE))
# nonsyn_total_1990 <- Reduce("+", nonsyn_list_1990) #95,168 variants
# nonsyn_list_2015 <- lapply(vcf_2015@region.data@synonymous, function(x) sum(x==0, na.rm=TRUE))
# nonsyn_total_2015 <- Reduce("+", nonsyn_list_2015) #84,644 variants

#Could also apply the same strategy for SNPs in various regions according to the GFF file
#Below are some examples for exon SNPs, but in principle could do the same thing for intron and UTR SNPs (if they're annotated in the GFF explicitly)
# exon_SNP_list_1990 <- lapply(vcf_1990@region.data@ExonSNPS, function(x) sum(x==1, na.rm=TRUE))
# exon_SNP_total_1990 <- Reduce("+", exon_SNP_list_1990) #153,487 variants


####Other misc notes####
#If I wanted to do analysis with the years combined, and the two populations as the years, I could do something like this to examine shared and fixed sites in the two years
#Read in the combined vcf, and then assign the two years as populations
#vcf_combined <- calc.fixed.shared(vcf_combined, fixed.threshold=.8) ?? could be something else
#vcf_combined@n.fixed.sites
#vcf_combined@n.shared.sites
#vcf_combined@n.monomorphic.sites


#######readVCF de-bugging notes#########

#When I tried to use readVCF to read in one contig at a time, I got some errors
#Neither of these worked, some error about a malformed GT field. Perhaps where data is missing?
#vcf_000000F <- readVCF("1990_gz/000000F.vcf.gz", numcols = 10000, tid = "000000F", frompos = 1, topos = 10000, approx = FALSE, gffpath = "GFF/000000F.gff")
#vcf_000001F <- readVCF("1990_gz/000001F.vcf.gz", numcols = 10000, tid = "000001F", frompos = 1, topos = 10000, approx = FALSE, gffpath = "GFF/000001F.gff")

#My thought was to use lapply or something to do this for all vcf files and make a big list of GENOME.class object, but I kept getting this error message for any contig I tested:
#GFF information ...
#vcff::open : file opened, contains 30 samples
#[1] "Available ContigIdentifiers (parameter tid):"
#[1] "000000F"
#VCF_readIntoCodeMatrix :: Malformed GT field!
#  Error in 1:numusedcols : NA/NaN argument
#In addition: Warning message:
#  In readVCF("1990_gz/000000F.vcf.gz", numcols = 10000, tid = "000000F",  :
#               NAs introduced by coercion

#So, I decided to try this with a short test file of 000000F.vcf 

#Lines with and without missing data for individuals (represented with .:.:.:.:.:.:.:.)
#vcf_test <- readVCF("read_vcf_test/test.vcf.gz", numcols = 10000, tid = "000000F", frompos = 1, topos = 10000, approx = FALSE, gffpath = "GFF/000000F.gff")

#Without any lines with missing data
#vcf_test2 <- readVCF("read_vcf_test/test2.vcf.gz", numcols = 10000, tid = "000000F", frompos = 1, topos = 10000, approx = FALSE, gffpath = "GFF/000000F.gff")

#Both work!! So what gives??

#If I use a small subset that has the missing data field at the END of a line, I get the same error as above
#vcf_test3 <- readVCF("read_vcf_test/test3.vcf.gz", numcols = 10000, tid = "000000F", frompos = 1, topos = 10000, approx = FALSE, gffpath = "GFF/000000F.gff")

#But if I take the test3 and remove only variant records that have missing individuals at the end, then it works!
#vcf_test4 <- readVCF("read_vcf_test/test4.vcf.gz", numcols = 10000, tid = "000000F", frompos = 1, topos = 10000, approx = FALSE, gffpath = "GFF/000000F.gff")

#So, readVCF will not work with VCF files that have lines with missing data at the end of a line, so this won't work for my analysis.
#Update 3/19/2018 - I heard from the author of PopGenome, and he mentioned that the GitHub version of PopGenome has this issue fixed.


####Testing out various LD analysis within PopGenome####
#Function to make dataframe for Wall's B
# wallB.df <- function(x, year, gene_type=NULL) {
#   if (is.null(gene_type))
#     df <- data.frame("Wall's B" = x@Wall.B, "Year" = year, "VCF" = x@region.names)
#   else
#     df <- data.frame("Wall's B" = x@Wall.B, "Year" = year, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", x@region.names, perl=TRUE), "Gene.Type" = gene_type)
#   df
# }
# 
# #Function to make dataframe for Kelly's ZnS
# kellyZ.df <- function(x, year, gene_type=NULL) {
#   if (is.null(gene_type))
#     df <- data.frame("Kelly's ZnS" = x@Kelly.Z_nS, "Year" = year, "VCF" = x@region.names)
#   else
#     df <- data.frame("Kelly's ZnS" = x@Kelly.Z_nS, "Year" = year, "Contig:Start:Stop" = sub("(.+F)\\.vcf : (\\d+) - (\\d+)", "\\1:\\2:\\3", x@region.names, perl=TRUE), "Gene.Type" = gene_type)
#   df
# }

#Linkage disequilibirum analysis#

# The B statistic is the proportion of pairs of adjacent segregating sites that are congruent, i.e., that have consistent genealogies

#The method considers each of the possible pairs of adjacent segregating sites. In some cases, the adjacent sites have the same pattern of variation among the individuals within the sample. In these cases, the two adjacent sites are in complete disequilibrium. Such pairs of sites are called congruent pairs. Wall's B is the proportion of pairs that are congruent. Wall's B can be considered a measure of LD with values approaching 1 indicating extensive congruence among adjacent segregating sites.

# vcf_1990_linkage <- linkage.stats(vcf_1990)
# vcf_2015_linkage <- linkage.stats(vcf_2015)

#Save and re-load objects for later analysis as linkage.stats takes a while to run
#save(vcf_1990_linkage, file="vcf_1990_linkage.Rdata")
#save(vcf_2015_linkage, file="vcf_2015_linkage.Rdata")
# load("vcf_1990_linkage.Rdata")
# load("vcf_2015_linkage.Rdata")


#Use WallB.df function to extract values
# vcf_1990_wallB <- wallB.df(vcf_1990_linkage, "1990")
# vcf_2015_wallB <- wallB.df(vcf_2015_linkage, "2015")
# 
# vcf_1990_kellyZ <- kellyZ.df(vcf_1990_linkage, "1990")
# vcf_2015_kellyZ <- kellyZ.df(vcf_2015_linkage, "2015")

#Make a dataframe for plotting
# wallB_df <- rbind(vcf_1990_wallB, vcf_2015_wallB)
# colnames(wallB_df)[[1]] <- "Wall's B" #Somehow I needed to rename this column again
# 
# kellyZ_df <- rbind(vcf_1990_kellyZ, vcf_2015_kellyZ)
# colnames(wallB_df)[[1]] <- "Kelly's ZnS" #Somehow I needed to rename this column again

#Get average across genome for 2 years
# wallB_1990 <- mean(vcf_1990_wallB$pop.1, na.rm = TRUE) #0.2020971
# wallB_2015 <- mean(vcf_2015_wallB$pop.1, na.rm = TRUE) #0.2716006
# 
# kellyZ_1990 <- mean(vcf_1990_kellyZ$pop.1, na.rm = TRUE) #0.1646942
# kellyZ_2015 <- mean(vcf_2015_kellyZ$pop.1, na.rm = TRUE) #0.1879721

#Test if significant difference between two years using the Mann-Whitney-Wilcoxon (Wilcoxon rank sum) test
# pv_wallB <- wilcox.test(vcf_1990_wallB$pop.1, vcf_2015_wallB$pop.1) #p < 2.2e-16 
# 
# pv_kellyZ <- wilcox.test(vcf_1990_kellyZ$pop.1, vcf_2015_kellyZ$pop.1) #p = 0.0001014
#For both metrics, LD is higher genome-wide for 2015 isolates compared to 1990 isolates

#Moreover, the module R2.stats is designed for fast compution of the correlation coefficient r2
#this doesn't work either on MSI w/ ff, or on my local machine without, so I think the function might be broken
#vcf_1990_r2 <- calc.R2(vcf_1990)

#4-gamete test, this does not run on my local machine or MSI. Just stalls out and never progresses.
# vcf_1990_4_gamete <- recomb.stats(vcf_1990)
# vcf_2015_4_gamete <- recomb.stats(vcf_2015)
