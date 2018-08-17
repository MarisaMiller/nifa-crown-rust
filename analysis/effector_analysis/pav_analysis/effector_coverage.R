#This script is for calculating coverage of effectors
#Change to the appropriate working directory where all coverage files are stored
setwd("")

#Load appropriate libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gridGraphics)

#Place all files with .cov extension in vector
#Make sure to use double backslash to escape the period in the filename
effector_cov_files <- list.files(pattern="^[19].*\\.effector\\.cov$", full.names=TRUE)
haus_effector_cov_files <- list.files(pattern="^[19].*\\.haus_effector\\.cov$", full.names=TRUE)

##Read in files and rename columns
#Rename columns vector
covcols <- c("contig", "NULL", "type", "start", "end", rep("NULL", 3), "gene_ID", "depth", "number_bases_at_depth", "gene_size", "percent_gene_at_depth")

#Make a function to read in files
loadCovFile <- function(x) {
  #read in a coverage file, extract the sample name from the file, and add it as a column
  #also add in a column if 1990 or 2015 isolate
  df <- read.delim(x, header=FALSE, col.names=covcols, colClasses = c("character", "NULL", "character", "numeric", "numeric", rep("NULL", 3), "character", rep("numeric", 4)))
  df$sample_name <- sub("(.+)\\..+\\.cov", "\\1", basename(x))
  df$sample_year <- sub("([12590]{2}).+\\.cov", "\\1", basename(x))
  df
}

#Then use lapply to apply this function to all files and load in as list
effector_cov <- lapply(effector_cov_files, loadCovFile)
haus_effector_cov <- lapply(haus_effector_cov_files, loadCovFile)

#Then, make this list into one long dataframe
effector_cov_df <- do.call(rbind, effector_cov)
haus_effector_cov_df <- do.call(rbind, haus_effector_cov)

#Remove "all" rows as they're not useful for this analysis
effector_cov_df <- subset(effector_cov_df, contig != "all")
haus_effector_cov_df <- subset(haus_effector_cov_df, contig != "all")

#Missing effectors (using whole gene body, from start to stop codon, not just CDS)
number_missing_effectors_by_sample <- effector_cov_df %>% group_by(.dots = c("sample_name", "gene_ID")) %>%
  filter(depth < 5) %>%
  summarize(fraction_gene_missing = sum(percent_gene_at_depth)) %>%
  filter(fraction_gene_missing > 0.75) %>%
  count(sample_name) %>%
  arrange(desc(n))

#Avg 7.72
#Min - max: 4 - 20

number_missing_haus_effectors_by_sample <- haus_effector_cov_df %>% group_by(.dots = c("sample_name", "gene_ID")) %>%
  filter(depth < 5) %>%
  summarize(fraction_gene_missing = sum(percent_gene_at_depth)) %>%
  filter(fraction_gene_missing > 0.75) %>%
  count(sample_name) %>%
  arrange(desc(n))

#Avg 2
#Min - max: 1 - 3

#Summarize results with a heatmap of fractions of effectors with a depth less than 5
#Which effectors are most commonly missing?
missing_effectors_df <- effector_cov_df %>% group_by(.dots = c("sample_name", "gene_ID")) %>%
  filter(sample_name != "12SD80") %>%
  filter(sample_name != "12NC29") %>%
  summarize(fraction_gene_missing = sum(if_else(depth < 5, percent_gene_at_depth, 0)))

missing_haus_effectors_df <- haus_effector_cov_df %>% group_by(.dots = c("sample_name", "gene_ID")) %>%
  filter(sample_name != "12SD80") %>%
  filter(sample_name != "12NC29") %>%
  summarize(fraction_gene_missing = sum(if_else(depth < 5, percent_gene_at_depth, 0)))

#Save tables for effector analysis
write.table(missing_effectors_df, file = "fraction_missing_effectors.txt", sep = "\t", row.names = FALSE)
write.table(missing_haus_effectors_df, file = "fraction_missing_haus_effectors.txt", sep = "\t", row.names = FALSE)

#Make matrix for heatmap, just for haustorial
#leave out the first two columns since they doesn't belong in the heatmap itself
wide_missing_haus_effectors_df <- spread(missing_haus_effectors_df, sample_name, fraction_gene_missing)

haus_effector_matrix <- as.matrix(wide_missing_haus_effectors_df[,c(2:ncol(wide_missing_haus_effectors_df))])

#save gene names for annotating the heatmap later if desired, and change to factor
row_names <- as.factor(wide_missing_haus_effectors_df$gene_ID)
rownames(haus_effector_matrix) <- row_names

#Reorder rows to be same as order in gff file
haus_gff <- read.delim("haus_expr_secreted_effectors_sd80_p.sorted.gff3", header = FALSE, colClasses=c("character", "NULL", "NULL", "numeric", "numeric", rep("NULL", 3), "character"), col.names = c("contig", "NULL", "NULL", "start", "end", rep("NULL", 3), "gene_ID"))

haus_effector_matrix <- haus_effector_matrix[match(haus_gff$gene_ID, rownames(haus_effector_matrix)),]

#change column names to be just the year and state rather than the isolate name
#colnames(haus_effector_matrix) <- sub("([15]{2})(\\w{2}).+", "20\\1_\\2", colnames(haus_effector_matrix), perl=TRUE)
#colnames(haus_effector_matrix) <- sub("([90]{2})(\\w{2}).+", "19\\1_\\2", colnames(haus_effector_matrix), perl=TRUE)

#Rearrange columns
#haus_effector_matrix <- haus_effector_matrix[, order(-haus_effector_matrix[which(rownames(haus_effector_matrix) == "ID=PCA_SD_18894"), ]) ]


#Now, make heatmap

#annotate columns by year
anno_df <- as.data.frame(ifelse(sub("([1590]{2}).+", "\\1", colnames(haus_effector_matrix)) == "15", "2015", "1990"))
colnames(anno_df) <- "Year"

ha <- HeatmapAnnotation(df = anno_df, col = list(Year = c("2015" =  "red", "1990" = "black")))

#Highlight the obviously deleted genes
subset <- which(rownames(haus_effector_matrix) %in% c("ID=PCA_SD_18894", "ID=PCA_SD_18895"))
labels <- rownames(haus_effector_matrix[88:89,])

#Different clustering methods with hclust
#The default method for clustering_distance_rows is "euclidean" and the default method for clustering_method_rows is "complete".
haus_plot <- Heatmap(haus_effector_matrix,
                        top_annotation = ha,
                        show_column_names = FALSE,
                        show_row_names = FALSE,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        width = unit(100, "mm"),
                        heatmap_legend_param = list(title = "Missing Gene Fraction", at = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), labels = c("0", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1")), col = colorRamp2(c(1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0), heat.colors(11))) + rowAnnotation(link = row_anno_link(at = subset, labels = labels), width = unit(.5, "cm") + max_text_width(labels))

#haus_plot <- draw(haus_plot, heatmap_legend_side = "left", annotation_legend_side = "top")

# Save the plots to a file
pdf("missing_haus_effectors.pdf", height = 8, width = 8)
haus_plot
dev.off()

#Save plot as R object to use in larger summary figure
save(haus_plot, file="haus_cov_plot.Rdata")

#Is there any obvious association with host status?

#annotate columns by year and host status
bt_data <- read.csv("../../buckthorn_analysis/buckthorn_FirstDateStateCountyCount.csv", header = TRUE)
no_bt_states <- bt_data %>% filter(PercentCounties < 10) %>% select(StateName)
bt_states <- bt_data %>% filter(PercentCounties > 10) %>% select(StateName)
mn_bt_isolates <- c("90MN1B-1", "90MN2B-1", "90MN3B-1", "90MN5B-1", "90MN7B-1", "90MN8B-2", "90MN9B-4", "90MN13B-3", "90MN14B-1", "90MN17B-1", "15MN13-4", "15MN14-4", "15MN15-3", "15MN16-3", "15MN17-5", "15MN18-1", "15MN18-3")

anno_host <- ifelse(sub("[1590]{2}(\\w{2}).+", "\\1", colnames(haus_effector_matrix), perl=TRUE) %in% no_bt_states$StateName, "NO BT", "BT")
anno_host[colnames(haus_effector_matrix) %in% mn_bt_isolates] <- "MN BT"
anno_host_df <- as.data.frame(anno_host)
colnames(anno_host_df) <- "Host Status"

ha_host <- HeatmapAnnotation(df = anno_host_df, col = list(`Host Status` = c("NO BT" =  "#E69F00", "BT" = "#D55E00", "MN BT" = "#56B4E9")))

haus_host_plot <- Heatmap(haus_effector_matrix,
                     name = "haus_host_plot",
                     top_annotation = ha_host,
                     bottom_annotation = ha,
                     show_column_names = FALSE,
                     show_row_names = FALSE,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     width = unit(100, "mm"),
                     heatmap_legend_param = list(title_position = "topcenter", title = "Missing Gene Fraction", at = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), labels = c("0", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1")), col = colorRamp2(c(1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0), heat.colors(11))) + rowAnnotation(link = row_anno_link(at = subset, labels = labels), width = unit(.5, "cm") + max_text_width(labels))

haus_host_plot <- draw(haus_host_plot, heatmap_legend_side = "left", annotation_legend_side = "top")

# decorate_heatmap_body("haus_host_plot", {
#   grid.lines(c(31/60, 31/60), c(0, 1), gp = gpar(col = "black", lwd = 2))
#   grid.lines(c(19/60, 19/60), c(0, 1), gp = gpar(col = "red", lwd = 1.5, lty = 2))
#   grid.lines(c(42/60, 42/60), c(0, 1), gp = gpar(col = "red", lwd = 1.5, lty = 2))
# })

# Save the plots to a file
pdf("missing_haus_effectors_withPop.pdf", height = 8, width = 8)
haus_host_plot
dev.off()

#What does the coverage look like in the 50kb windows around this gene?
#Import file from the genome-wide window coverage analysis
window_cov <- read.delim("../../coverage_analysis/genome_window_analysis/fraction_missing_per_window.txt", header=TRUE)

window_cov_000335_subset <- window_cov %>% filter(contig_window == "000335F_1" & start == 0 & end == 50000)

#What's the average coverage across the window in the two isolate groups either lacking or containing SD_18894/95?
isolates_containing_18894_5 <- missing_haus_effectors_df %>% filter(gene_ID == "ID=PCA_SD_18894" | gene_ID == "ID=PCA_SD_18895", fraction_gene_missing == 0) %>% select(sample_name)

isolates_lacking_18894_5 <- missing_haus_effectors_df %>% filter(gene_ID == "ID=PCA_SD_18894" | gene_ID == "ID=PCA_SD_18895", fraction_gene_missing > 0) %>% select(sample_name)

mean_window_cov_with18894_5 <- window_cov_000335_subset %>% filter(sample_name %in% isolates_containing_18894_5$sample_name) %>% summarize(mean(fraction_window_missing), sd(fraction_window_missing))
# mean(fraction_window_missing) sd(fraction_window_missing)
#                     0.01266069                 0.006723293
#Missing ~630 bp

mean_window_cov_without18894_5 <- window_cov_000335_subset %>% filter(sample_name %in% isolates_lacking_18894_5$sample_name) %>% summarize(mean(fraction_window_missing), sd(fraction_window_missing))
# mean(fraction_window_missing) sd(fraction_window_missing)
#                      0.1171935                  0.00632534
#Missing ~5.8 kb
