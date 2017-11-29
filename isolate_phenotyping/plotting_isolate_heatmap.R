#ComplexHeatmap_1.14.0 and colorRamps_2.3 used to make heatmaps

#set working directory to where files are
setwd("")
list.files()

#load required packages
library(ComplexHeatmap)
library(circlize)
#library(RColorBrewer) #another alternative package for colors if so desired
library(colorRamps)

#read in data and inspect
scoring_data = read.csv("2015_phenotypes.csv", header=TRUE)
#scoring_data = read.csv("1990_phenotypes.csv", header=TRUE)

names(scoring_data) <- c("Differential Line", "15ND19-2",	"15ND19-5",	"15ND20-3",	"15ND20-4",	"15MN13-4",	"15MN14-4",	"15MN15-3",	"15MN16-3",	"15MN17-5",	"15MN18-1",	"15MN18-3",	"15MN10-4",	"15MN10-5",	"15MN27-3",	"15MN23-1",	"15MN24-1",	"15MN25-3",	"15SD30-1",	"15SD30-3",	"15SD11-1",	"15SD11-2",	"15NE8-4",	"15NE8-5",	"15NE9-1",	"15NE9-3",	"15OH12-3",	"15MS7-1",	"15TX3-1",	"15FL1-2",	"15FL1-4")
#names(scoring_data) <- c("Differential Line", "90MN1B-1",	"90MN2B-1",	"90MN3B-1",	"90MN5B-1",	"90MN7B-1",	"90MN8B-2",	"90MN9B-4",	"90MN13B-3",	"90MN14B-1",	"90MN17B-1",	"90MN137-1",	"90MN148-1",	"90MN149-2",	"90MN152-1",	"90MN153-1",	"90WI131-1",	"90WI132-1",	"90SD164-1",	"90SD171-1",	"90SD172-1",	"90PA162-1",	"90KS101-1",	"90AR100-1",	"90GA16-1",	"90LA38-1",	"90TX45-1",	"90TX47-1",	"90TX52-1",	"90TX58-1",	"90TX70-1")

head(scoring_data)

#make the heatmap data into a matrix
#Use absolute values because of large number of annotations
scoring_matrix <- as.matrix(scoring_data[  ,c(2:ncol(scoring_data))])

#[all rows, columns x-x]
#leave out the first column since it doesn't belong in the heatmap itself

#save TF names for annotating the heatmap later
row_names <- scoring_data$"Differential Line"
head(row_names)
class(row_names) # Needs to be factor and not dataframe to name matrix rows

rownames(scoring_matrix) = row_names

#Now, make heatmap

#Different clustering methods with hclust
#The default method for clustering_distance_rows is "euclidean" and the default method for clustering_method_rows is "complete".
scoring_plot <- Heatmap(scoring_matrix,
        row_names_side = "left",
        column_names_side = "top",
        row_names_gp = gpar(cex=0.8, fontface = ifelse(grepl("Pc", row_names), "italic", "plain")),
        column_names_gp = gpar(cex=0.8),
        width = unit(100, "mm"),
        heatmap_legend_param = list(title = NULL, at = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")), col = colorRamp2(c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0), heat.colors(10)))
#col = colorRamp2(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), blue2green(10))) #another color scheme

# Save the plots to a file
pdf("scoring_plot_2015.pdf")
#pdf("scoring_plot_1990.pdf")
scoring_plot
dev.off()
