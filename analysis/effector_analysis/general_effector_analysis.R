#This script is for further analysis of effectors and to summarize variants, variant annotations, and other population genetic summary statistics.

#Change working directory to appropriate location
setwd(")"

#Load appropriate libraries
library(GenomicFeatures)
library(rtracklayer)
library(ggbio)
library(dplyr) #Make sure to load dplyr after GenomicFeatures to avoid issues with the select function
library(forcats)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

#First read in coverage data to begin building summary table
pav_all_effectors <- read.delim("./pav_analysis/fraction_missing_effectors.txt", header = TRUE)

#Summarize coverage data by classifying the percent of isolates with >0.75 of the gene missing
wide_pav_all_effectors <- spread(pav_all_effectors, sample_name, fraction_gene_missing)
wide_pav_all_effectors$fraction_isolates_missing <- (rowSums(wide_pav_all_effectors[,2:ncol(wide_pav_all_effectors)] > 0.75)/60)

effector_summary <- wide_pav_all_effectors[,c(1,62)]

#Add population genetic data to summary table and add all to larger summary table
nucdiv_all_effectors <- read.delim("../popgenome/nucdiv_effectors.txt", header = TRUE)
nucdiv_all_effectors <- nucdiv_all_effectors %>% dplyr::select(-Contig.Start.Stop)
effector_summary <- left_join(effector_summary, nucdiv_all_effectors, by = c("gene_ID"))

theta_all_effectors <- read.delim("../popgenome/theta_effectors.txt", header = TRUE)
theta_all_effectors <- theta_all_effectors %>% dplyr::select(-Contig.Start.Stop)
effector_summary <- left_join(effector_summary, theta_all_effectors, by = c("gene_ID"))

TajD_all_effectors <- read.delim("../popgenome/TajD_effectors.txt", header = TRUE)
TajD_all_effectors <- TajD_all_effectors %>% dplyr::select(-Contig.Start.Stop)
effector_summary <- left_join(effector_summary, TajD_all_effectors, by = c("gene_ID"))

Fst_all_effectors <- read.delim("../popgenome/Fst_combined_effectors.txt", header = TRUE)
Fst_all_effectors <- Fst_all_effectors %>% dplyr::select(-c(Contig.Start.Stop, Year))
effector_summary <- left_join(effector_summary, Fst_all_effectors, by = c("gene_ID"))

clr_all_effectors <- read.delim("../popgenome/clr_effectors.txt", header = TRUE)
clr_all_effectors <- clr_all_effectors %>% dplyr::select(-c(Contig.Start.Stop))
effector_summary <- left_join(effector_summary, clr_all_effectors, by = c("gene_ID"))

#Clean up column names and do a little formatting
effector_summary$fraction_isolates_missing <- round(effector_summary$fraction_isolates_missing, 2)
effector_summary[,2:ncol(effector_summary)] <- round(effector_summary[,2:ncol(effector_summary)], 5)

colnames(effector_summary) <- c("Gene ID", "Fraction of isolates with >0.75 of gene missing", "Nucleotide Diversity 1990", "Nucleotide Diversity 2015", "Watterson's Theta 1990", "Watterson's Theta 2015", "Tajima's D 1990", "Tajima's D 2015", "Fst 1990 vs 2015", "Composite Likelihood Ratio 1990", "Composite Likelihood Ratio 2015")

#Save to table for publication
write.table(effector_summary, file = "effector_summary.txt", sep = "\t", row.names = FALSE)


#Make summary plots for various statistics

#make a new summary to include all isolate data for coverage (rather than summarized for all isolates for each gene)
effector_summary_plot_format <- effector_summary %>% dplyr::select(-`Fraction of isolates with >0.75 of gene missing`) %>% left_join(., rename(wide_pav_all_effectors, "Gene ID" = gene_ID), by = "Gene ID") %>% dplyr::select(-fraction_isolates_missing)

#Subset to just haustorial effectors for plot
haus_gff <- read.delim("./pav_analysis/haus_expr_secreted_effectors_sd80_p.sorted.gff3", header = FALSE, colClasses=c("character", "NULL", "NULL", "numeric", "numeric", rep("NULL", 3), "character"), col.names = c("contig", "NULL", "NULL", "start", "end", rep("NULL", 3), "Gene ID"))

effector_summary_plot_format_haus <- left_join(rename(haus_gff, "Gene ID" = Gene.ID), effector_summary_plot_format, by = "Gene ID")

#Removing the line below will leave the NAs in the table rather than change to 0
#This is preferable because some statistics (like nucleotide diversity) can actually be 0, but some stats (like Tajima's D) can be NA if they can't be calculated due to 0 nucleotide diversity in a given year (or both years) for a given gene
#effector_summary_plot_format_haus[, 5:13][is.na(effector_summary_plot_format_haus[, 5:13])] <- 0

#Break into 3 separate dfs for easier plotting
nuc_div_haus <- effector_summary_plot_format_haus[,c(4,5,6)]
nuc_div_haus$change <- abs(nuc_div_haus$`Nucleotide Diversity 2015` - nuc_div_haus$`Nucleotide Diversity 1990`) #Add a column of the abs value of the difference between the two years
nuc_div_haus$`Gene ID` <- factor(nuc_div_haus$`Gene ID`)
nuc_div_haus <- nuc_div_haus %>% arrange(desc(change))
#Rearrange levels of factor based on change
nuc_div_haus$`Gene ID` <- fct_rev(factor(nuc_div_haus$`Gene ID`, levels=nuc_div_haus$`Gene ID`[order(nuc_div_haus$change)], ordered=TRUE))
nuc_div_haus_long <- nuc_div_haus %>% dplyr::select(-change) %>% gather(variable, "Nucleotide Diversity", -`Gene ID`)
nuc_div_haus_long$year <- sub(".+([19025]{4})", "\\1", nuc_div_haus_long$variable)
nuc_div_haus_long$variable <- sub("(.+)[19025]{4}", "\\1", nuc_div_haus_long$variable)

#Order Tajima's D, CLR, and Fst based on the magnitutde of change in the nucleotide diversity to keep gene order consistent so plots can be directly compared 
tajd_haus <- effector_summary_plot_format_haus[,c(4,9,10)]
tajd_haus$change <- abs(tajd_haus$`Tajima's D 2015` - tajd_haus$`Tajima's D 1990`) #Add a column of the abs value of the difference between the two years
tajd_haus$`Gene ID` <- factor(tajd_haus$`Gene ID`)
tajd_haus <- tajd_haus %>% arrange(desc(change))
#Rearrange levels of factor based on change in nuc div
tajd_haus$`Gene ID` <- fct_rev(factor(tajd_haus$`Gene ID`, levels=nuc_div_haus$`Gene ID`[order(nuc_div_haus$change)], ordered=TRUE))
tajd_haus_long <- tajd_haus %>% dplyr::select(-change) %>% gather(variable, "Tajima's D", -`Gene ID`)
tajd_haus_long$year <- sub(".+([19025]{4})", "\\1", tajd_haus_long$variable)
tajd_haus_long$variable <- sub("(.+)[19025]{4}", "\\1", tajd_haus_long$variable)

clr_haus <- effector_summary_plot_format_haus[,c(4,12,13)]
clr_haus$change <- abs(clr_haus$`Composite Likelihood Ratio 2015` - clr_haus$`Composite Likelihood Ratio 1990`) #Add a column of the abs value of the difference between the two years
clr_haus$`Gene ID` <- factor(clr_haus$`Gene ID`)
clr_haus <- clr_haus %>% arrange(desc(change))
#Rearrange levels of factor based on change in nuc div
clr_haus$`Gene ID` <- fct_rev(factor(clr_haus$`Gene ID`, levels=nuc_div_haus$`Gene ID`[order(nuc_div_haus$change)], ordered=TRUE))
clr_haus_long <- clr_haus %>% dplyr::select(-change) %>% gather(variable, "Composite Likelihood Ratio", -`Gene ID`)
clr_haus_long$year <- sub(".+([19025]{4})", "\\1", clr_haus_long$variable)
clr_haus_long$variable <- sub("(.+)[19025]{4}", "\\1", clr_haus_long$variable)

fst_haus <- effector_summary_plot_format_haus[,c(4,11)]
fst_haus$`Gene ID` <- factor(fst_haus$`Gene ID`)
fst_haus <- fst_haus %>% arrange(desc(`Fst 1990 vs 2015`))
#Rearrange levels of factor based on change in nuc div
fst_haus$`Gene ID` <- fct_rev(factor(fst_haus$`Gene ID`, levels=nuc_div_haus$`Gene ID`[order(nuc_div_haus$change)], ordered=TRUE))

#Set theme to be classic
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

nuc_div_plot <- ggplot(nuc_div_haus_long, aes(x = `Gene ID`, y = `Nucleotide Diversity`)) +
  geom_line(color = "grey") +
  geom_point(aes(color = nuc_div_haus_long$year)) +
  scale_color_manual(values=c("black", "red"), guide=FALSE) +
  theme(axis.text.x = element_blank()) +
  xlab("")

tajd_plot <- ggplot(tajd_haus_long, aes(x = `Gene ID`, y = `Tajima's D`)) +
  geom_line(color = "grey") +
  geom_point(aes(color = tajd_haus_long$year)) +
  scale_color_manual(values=c("black", "red"), guide=FALSE) +
  theme(axis.text.x = element_blank()) +
  xlab("")

#See which effectors are above 95th percentile for the CLR value and add lines to plot
quantile(effector_summary$`Composite Likelihood Ratio 1990`, 0.95, na.rm = TRUE)
# 95% 
# 82.8901 

quantile(effector_summary$`Composite Likelihood Ratio 2015`, 0.95, na.rm = TRUE)
# 95% 
# 86.80082 

clr_plot <- ggplot(clr_haus_long, aes(x = `Gene ID`, y = `Composite Likelihood Ratio`)) +
  geom_line(color = "grey") +
  geom_point(aes(color = tajd_haus_long$year)) +
  scale_color_manual(values=c("black", "red"), guide=FALSE) +
  theme(axis.text.x = element_blank()) +
  xlab("") +
  geom_hline(yintercept=82.8901, color = "black", size=0.5) +
  geom_hline(yintercept=86.80082, color = "red", size=0.5)

fst_plot <- ggplot(fst_haus, aes(x = `Gene ID`, y = `Fst 1990 vs 2015`)) +
  geom_line(color = "grey") +
  geom_point() +
  theme(axis.text.x = element_blank()) +
  xlab("Haustorial Effectors") +
  ylab(expression(italic("F")[ST]))


#Need to account for differing numbers of y-axis characters shifting plots out of alignment
gA <- ggplotGrob(nuc_div_plot)
gB <- ggplotGrob(tajd_plot)
gC <- ggplotGrob(clr_plot)
gD <- ggplotGrob(fst_plot)
maxWidth = grid::unit.pmax(gA$widths[3:5], gB$widths[3:5], gC$widths[3:5], gD$widths[3:5])
gA$widths[3:5] <- as.list(maxWidth)
gB$widths[3:5] <- as.list(maxWidth)
gC$widths[3:5] <- as.list(maxWidth)
gD$widths[3:5] <- as.list(maxWidth)

pdf("haus_effector_summary.pdf", width = 7, height = 9)
grid.arrange(gA, gB, gC, gD, nrow = 4)
dev.off()

#Coverage plot
cov_effector_matrix <- as.matrix(effector_summary_plot_format_haus[,c(which(names(effector_summary_plot_format_haus)=="15FL1-2"):ncol(effector_summary_plot_format_haus))])

#save gene names for annotating the heatmap later if desired
row_names <- as.factor(effector_summary_plot_format_haus$`Gene ID`)
rownames(cov_effector_matrix) <- row_names

#Reorder to have 1990 isolates first and sort by row means
cov_effector_matrix <- cov_effector_matrix[,c(31:ncol(cov_effector_matrix), 1:30)]
cov_effector_matrix <- cov_effector_matrix[order(rowMeans(cov_effector_matrix), decreasing = T), ]

#Annotate by year
cov_anno_df <- as.data.frame(ifelse(sub("([1590]{2}).+", "\\1", colnames(cov_effector_matrix)) == "15", "2015", "1990"))
colnames(cov_anno_df) <- "Year"
cov_ha <- HeatmapAnnotation(df = cov_anno_df, col = list(Year = c("2015" =  "red", "1990" = "black")))

#Highlight the obviously deleted genes
cov_subset <- which(rownames(cov_effector_matrix) %in% c("ID=PCA_SD_18894", "ID=PCA_SD_18895"))
cov_labels <- rownames(cov_effector_matrix[1:2,])

#Now, make heatmap
cov_haus_plot <- Heatmap(cov_effector_matrix,
                         top_annotation = cov_ha,
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         width = unit(80, "mm"),
                         rect_gp = gpar(col = "lightgrey", lwd = 0.05),
                         heatmap_legend_param = list(title = "Missing Gene Fraction", at = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), labels = c("0", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1")), col = colorRamp2(c(1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0), heat.colors(11))) + rowAnnotation(link = row_anno_link(at = cov_subset, labels = cov_labels), width = unit(.5, "cm") + max_text_width(cov_labels))

# Save the plot to a file
pdf("coverage_haus_plot.pdf", height = 7, width = 10)
cov_haus_plot 
dev.off()


#########Make ggbio coverage summary across region with SD_18894 and SD_18895##########

#Use chrominfo files when making TxDb. Not sure if needed but doesn't hurt to include
chromlen <- read.delim("Puccinia_coronata_avenae_12SD80.primary.chromlen", header = FALSE)
colnames(chromlen) <- c("chrom", "length")
chromlen$is_circular <- FALSE


#There was an issue with the stop and start coordinates for genes on the negative strand, so I removed those offending lines and subsetted the gff3 to just the contig of interest for the plot
txdb <- makeTxDbFromGFF(file="Puccinia_coronata_avenae_12SD80.000335F.gff3",
                        chrominfo=chromlen,
                        dataSource="gff3 from 12SD80 000335F primary contig",
                        organism="Puccinia coronata")


#Read in repeats with rtracklayer
repeats_gr <- import("Puccinia_coronata_avenae_12SD80.000335F.repeats.gff3")

#Point just to chrom of interest
seqlevels(txdb) <- "000335F"

#Subset to just genes
wh <- subset(genes(txdb))

#Make plot of just gene track in the region of interest on 000335F
p.txdb <- ggbio::autoplot(txdb, which = wh) + xlim(30000, 50000)
repeats <- ggbio::autoplot(repeats_gr, which = wh) + xlim(30000, 50000)

#List all bam files
bam_list <- setNames(as.list(list.files(path = "./subset_alignments", pattern = "*.bam$", full.names = TRUE)), sprintf("isolate_%s",seq(1:60)))

#Make list of plot objects
bam_plots <- lapply(bam_list, function(x) ggbio::autoplot(x, which = wh) + xlim(30000, 50000))

#Give them easy names (isolate_1-60)
for(i in 1:length(bam_plots)) assign(names(bam_plots)[i], bam_plots[[i]])

#Code to make very large plot, although this is unreasonable for a figure
#tks <- tracks(isolate_1, isolate_2, isolate_3, isolate_4, isolate_5, isolate_6, isolate_7, isolate_8, isolate_9, isolate_10, isolate_11, isolate_12, isolate_13, isolate_14, isolate_15, isolate_16, isolate_17, isolate_18, isolate_19, isolate_20, isolate_21, isolate_22, isolate_23, isolate_24, isolate_25, isolate_26, isolate_27, isolate_28, isolate_29, isolate_30, isolate_31, isolate_32, isolate_33, isolate_34, isolate_35, isolate_36, isolate_37, isolate_38, isolate_39, isolate_40, isolate_41, isolate_42, isolate_43, isolate_44, isolate_45, isolate_46, isolate_47, isolate_48, isolate_49, isolate_50, isolate_51, isolate_52, isolate_53, isolate_54, isolate_55, isolate_56, isolate_57, isolate_58, isolate_59, isolate_60, gene = p.txdb, repeats)

# pdf("ggbio_test.pdf", height = 40)
# tks
# dev.off()


#Instead, just show two isolates from each class (w/ and w/o deletion) to illustrate what the deletion looks like
tks <- tracks(isolate_1, isolate_4, isolate_5, isolate_2, isolate_3, isolate_10, gene = p.txdb, repeats) + theme(axis.text=element_text(size=10))

pdf("deletion_region.pdf", height = 9)
tks
dev.off()

#I wanted to get this to work, but the text commands I built were just not interpreted correctly
# bam_command_list <- lapply(bam_list, function(x) noquote(paste("(ggbio::autoplot(", dQuote(x), ", which = wh) + xlim(30000, 40000))", sep="")))
#Make list to comma separated string. If I add unquote again here, then I get an error that tracks doesn't accept an unoquote object. If I leave it off, then there is an error about not finding file extensions, probably because the commands just aren't interpreted correctly.
# bam_command_string <- toString(bam_command_list, sep = ", ")
#The command below does work though, where I just copy and paste some of the output of the bam_string_command into tracks, so I know that this approach is almost correct.
# tracks((ggbio::autoplot('./subset_alignments/15FL1-2.000335F.bam', which = wh) + xlim(30000, 40000)), (ggbio::autoplot('./subset_alignments/15Fl1-4.000335F.bam', which = wh) + xlim(30000, 40000)), (ggbio::autoplot('./subset_alignments/15MN10-4.000335F.bam', which = wh) + xlim(30000, 40000)), (ggbio::autoplot('./subset_alignments/15MN10-5.000335F.bam', which = wh) + xlim(30000, 40000)), gene = p.txdb)


#Summarize exonic SNP and Indel variants per effector broken down by functional type.
#I chose to focus on these since they have the largest chance to impact gene function.

#Read in files generated using bedtools to intersect the Annovar output with a gff3 file containing coordinates of effectors
exonicEffectorVariants_1990 = read.delim("./variation_analysis/1990_annovar_exonic_effector_variants.txt", header=FALSE, colClasses = c("character", "numeric", "NULL", "character", "character", "NULL", "NULL", "character", "character", rep("character", 30), "NULL", "character", "character", rep("NULL", 5), "numeric", rep("NULL", 5), "numeric", "numeric", "NULL", "character", "NULL", "character"))

colnames(exonicEffectorVariants_1990) = c("contig", "position", "ref", "alt", "info", "format", "90GA16-1", "90TX52-1", "90TX47-1", "90MN2B-1", "90MN17B-1", "90MN149-2", "90MN153-1", "90MN14B-1", "90KS101-1", "90MN148-1", "90WI132-1", "90MN7B-1", "90MN13B-3", "90TX45-1", "90SD172-1", "90MN8B-2", "90WI131-1", "90MN9B-4", "90MN3B-1", "90SD171-1", "90MN152-1", "90MN1B-1", "90LA38-1", "90MN137-1", "90PA162-1", "90TX70-1", "90TX58-1", "90SD164-1", "90AR100-1", "90MN5B-1", "type", "gene:transcript_identifier:transcript_sequence_change", "allele_frequency", "start", "stop", "strand", "Gene ID")

exonicEffectorVariants_2015 = read.delim("./variation_analysis/2015_annovar_exonic_effector_variants.txt", header=FALSE, colClasses = c("character", "numeric", "NULL", "character", "character", "NULL", "NULL", "character", "character", rep("character", 25), "NULL", rep("character", 2), "NULL", rep("character", 3), "NULL", "character", "character", rep("NULL", 5), "numeric", rep("NULL", 5), "numeric", "numeric", "NULL", "character", "NULL", "character"))

colnames(exonicEffectorVariants_2015) = c("contig", "position", "ref", "alt", "info", "format", "15ND20-3", "15MN25-3", "15NE9-1", "15ND20-4", "15MN24-1", "15NE8-4", "15FL1-4", "15MN17-5", "15SD30-3", "15MN16-3", "15SD30-1", "15MN18-1", "15SD11-2", "15NE9-3", "15ND19-5", "15MN27-3", "15OH12-3", "15TX3-1", "15SD11-1", "15FL1-2", "15NE8-5", "15MN14-4", "15MN23-1", "15MN15-3", "15MN10-5", "15MN18-3", "15ND19-2", "15MS7-1", "15MN10-4", "15MN13-4", "type", "gene:transcript_identifier:transcript_sequence_change", "allele_frequency", "start", "stop", "strand", "Gene ID")

#Replace genotype information with 0 or 1 to mark if gene has variation in a given isolate
exonicEffectorVariants_1990[,7:36] <- ifelse(sapply(exonicEffectorVariants_1990[,7:36], grepl, pattern = "(\\d)/(1):.+", perl = TRUE), 1, 0)
exonicEffectorVariants_2015[,7:36] <- ifelse(sapply(exonicEffectorVariants_2015[,7:36], grepl, pattern = "(\\d)/(1):.+", perl = TRUE), 1, 0)

#Note that the counts are not necessarily accurate, but are a rough approximation not taking into account zygosity
#If a gene of interest is noted to have some variants, then we can go back to the original output from Annovar and see what is present.
summaryExonicEffectorVariants_1990 <- exonicEffectorVariants_1990 %>% gather(., key = Isolate, value = Variation, "90GA16-1", "90TX52-1", "90TX47-1", "90MN2B-1", "90MN17B-1", "90MN149-2", "90MN153-1", "90MN14B-1", "90KS101-1", "90MN148-1", "90WI132-1", "90MN7B-1", "90MN13B-3", "90TX45-1", "90SD172-1", "90MN8B-2", "90WI131-1", "90MN9B-4", "90MN3B-1", "90SD171-1", "90MN152-1", "90MN1B-1", "90LA38-1", "90MN137-1", "90PA162-1", "90TX70-1", "90TX58-1", "90SD164-1", "90AR100-1", "90MN5B-1") %>% group_by(.dots = c(as.name("Gene ID"), "Isolate", "type")) %>% summarize(Type_Count = sum(Variation)) %>% spread(Isolate, Type_Count)

summaryExonicEffectorVariants_2015 <- exonicEffectorVariants_2015 %>% gather(., key = Isolate, value = Variation, "15ND20-3", "15MN25-3", "15NE9-1", "15ND20-4", "15MN24-1", "15NE8-4", "15FL1-4", "15MN17-5", "15SD30-3", "15MN16-3", "15SD30-1", "15MN18-1", "15SD11-2", "15NE9-3", "15ND19-5", "15MN27-3", "15OH12-3", "15TX3-1", "15SD11-1", "15FL1-2", "15NE8-5", "15MN14-4", "15MN23-1", "15MN15-3", "15MN10-5", "15MN18-3", "15ND19-2", "15MS7-1", "15MN10-4", "15MN13-4") %>% group_by(.dots = c(as.name("Gene ID"), "Isolate", "type")) %>% summarize(Type_Count = sum(Variation)) %>% spread(Isolate, Type_Count)

#An example of how this could be added to larger summary table
#effector_summary <- left_join(effector_summary, summaryExonicEffectorVariants_1990, by = c("gene_ID"))


















#For heatmap versions of statistics - not nearly as effective
# #Sort by Tajima's D
# effector_summary_plot_format_haus <- effector_summary_plot_format_haus %>% arrange(desc(`Tajima's D 1990`), `Tajima's D 2015`)
# 
# #Variant heatmap
# variant_effector_matrix <- as.matrix(effector_summary_plot_format_haus[,5:14])
# 
# #Re-order columns to be more logical
# variant_effector_matrix <- variant_effector_matrix[,c(1,2,3,4,5,6,7,10,8,9)]
# 
# #save gene names for annotating the heatmap later if desired
# rownames(variant_effector_matrix) <- row_names
# 
# #Now, make heatmap
# variant_haus_plot <- Heatmap(variant_effector_matrix,
#                      show_column_names = TRUE,
#                      column_names_gp = gpar(fontsize = 8),
#                      column_names_side = "bottom",
#                      show_row_names = FALSE,
#                      cluster_rows = FALSE,
#                      cluster_columns = FALSE,
#                      width = unit(30, "mm"),
#                      col = colorRamp2(c(0, 20), c("lightgoldenrodyellow", "black")),
#                      heatmap_legend_param = list(title = "Variant Count"))
# 
# #Nucleotide diversity heatmap
# nucdiv_effector_matrix <- as.matrix(effector_summary_plot_format_haus[,15:16])
# 
# #save gene names for annotating the heatmap later if desired
# rownames(nucdiv_effector_matrix) <- row_names
# 
# #Annotate by year
# nucdiv_anno_df <- as.data.frame(ifelse(sub(".+([21590]{4})", "\\1", colnames(nucdiv_effector_matrix)) == "2015", "2015", "1990"))
# colnames(nucdiv_anno_df) <- "Year"
# 
# nucdiv_ha <- HeatmapAnnotation(df = nucdiv_anno_df, col = list(Year = c("2015" =  "red", "1990" = "black")))
# 
# #Now, make heatmap
# nucdiv_haus_plot <- Heatmap(nucdiv_effector_matrix,
#                              top_annotation = nucdiv_ha,
#                              show_column_names = FALSE,
#                              show_row_names = FALSE,
#                              cluster_rows = FALSE,
#                              cluster_columns = FALSE,
#                              width = unit(20, "mm"),
#                              col = colorRamp2(c(0, .015), c("lightgoldenrodyellow", "blue")),
#                              heatmap_legend_param = list(title = "Nucleotide Diversity"))
# 
# #Tajima's D heatmap
# tajD_effector_matrix <- as.matrix(effector_summary_plot_format_haus[,19:20])
# 
# #save gene names for annotating the heatmap later if desired
# rownames(tajD_effector_matrix) <- row_names
# 
# #Annotate by year
# tajD_anno_df <- as.data.frame(ifelse(sub(".+([21590]{4})", "\\1", colnames(tajD_effector_matrix)) == "2015", "2015", "1990"))
# colnames(tajD_anno_df) <- "Year"
# 
# tajD_ha <- HeatmapAnnotation(df = tajD_anno_df, col = list(Year = c("2015" =  "red", "1990" = "black")))
# 
# #Now, make heatmap
# tajD_haus_plot <- Heatmap(tajD_effector_matrix,
#                             top_annotation = tajD_ha,
#                             show_column_names = FALSE,
#                             show_row_names = FALSE,
#                             cluster_rows = FALSE,
#                             cluster_columns = FALSE,
#                             width = unit(20, "mm"),
#                             col = colorRamp2(c(-3, 0, 4), c("orange3", "white", "purple4")),
#                             heatmap_legend_param = list(title = "Tajima's D"))
# 
# #Fst heatmap
# fst_effector_matrix <- as.matrix(effector_summary_plot_format_haus[,21])
# 
# #save gene names for annotating the heatmap later if desired
# rownames(fst_effector_matrix) <- row_names
# 
# #Now, make heatmap
# fst_haus_plot <- Heatmap(fst_effector_matrix,
#                           show_column_names = FALSE,
#                           show_row_names = FALSE,
#                           cluster_rows = FALSE,
#                           cluster_columns = FALSE,
#                           width = unit(15, "mm"),
#                           col = colorRamp2(c(-0.03, 0, .5), c("saddlebrown", "white", "seagreen4")),
#                           heatmap_legend_param = list(title = expression("F"[st])))
# 
# 
# ht_list = variant_haus_plot + fst_haus_plot + tajD_haus_plot + nucdiv_haus_plot + cov_haus_plot
# 
# haus_plot <- draw(ht_list, heatmap_legend_side = "left", annotation_legend_side = "left")
# 
# # Save the plot to a file
# pdf("summary_haus_plot.pdf", height = 9, width = 12)
# haus_plot 
# dev.off()
# 
# 
# #Version without variant annotation
# ht_list = fst_haus_plot + tajD_haus_plot + nucdiv_haus_plot + cov_haus_plot
# 
# haus_plot <- draw(ht_list, heatmap_legend_side = "left", annotation_legend_side = "left")
# 
# # Save the plot to a file
# pdf("summary_haus_plot_noAnnovar.pdf", height = 9, width = 12)
# haus_plot 
# dev.off()
