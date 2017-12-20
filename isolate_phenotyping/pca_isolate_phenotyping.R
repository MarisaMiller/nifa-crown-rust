#Use pcaMethods to conduct PCA analysis of isolate phenotyping data
#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaMethods")
library(pcaMethods)

#set working directory to where files are, modify the path
setwd("")

#read in data and inspect
scoring_data_2015 = read.csv("2015/2015_phenotypes.csv", header=FALSE)
scoring_data_1990 = read.csv("1990/1990_phenotypes.csv", header=FALSE)

#Remove first row
scoring_data_2015 = scoring_data_2015[-1,]
scoring_data_1990 = scoring_data_1990[-1,]

#Save differential names to add back after transposition
differentials_2015 <- scoring_data_2015$V1
differentials_1990 <- scoring_data_1990$V1

#Transpose, removing first column while doing so
transposed_scoring_data_2015 <- as.data.frame(t(scoring_data_2015[,-1]))
transposed_scoring_data_1990 <- as.data.frame(t(scoring_data_1990[,-1]))

#Add back in columnn names
colnames(transposed_scoring_data_2015) <- differentials_2015
colnames(transposed_scoring_data_1990) <- differentials_1990

#Add in rownames
rownames(transposed_scoring_data_2015) <- c("15ND19-2",	"15ND19-5",	"15ND20-3",	"15ND20-4",	"15MN13-4",	"15MN14-4",	"15MN15-3",	"15MN16-3",	"15MN17-5",	"15MN18-1",	"15MN18-3",	"15MN10-4",	"15MN10-5",	"15MN27-3",	"15MN23-1",	"15MN24-1",	"15MN25-3",	"15SD30-1",	"15SD30-3",	"15SD11-1",	"15SD11-2",	"15NE8-4",	"15NE8-5",	"15NE9-1",	"15NE9-3",	"15OH12-3",	"15MS7-1",	"15TX3-1",	"15FL1-2",	"15FL1-4")
rownames(transposed_scoring_data_1990) <- c("90MN1B-1",	"90MN2B-1",	"90MN3B-1",	"90MN5B-1",	"90MN7B-1",	"90MN8B-2",	"90MN9B-4",	"90MN13B-3",	"90MN14B-1",	"90MN17B-1",	"90MN137-1",	"90MN148-1",	"90MN149-2",	"90MN152-1",	"90MN153-1",	"90WI131-1",	"90WI132-1",	"90SD164-1",	"90SD171-1",	"90SD172-1",	"90PA162-1",	"90KS101-1",	"90AR100-1",	"90GA16-1",	"90LA38-1",	"90TX45-1",	"90TX47-1",	"90TX52-1",	"90TX58-1",	"90TX70-1")

#Let's simplify for the sake of the plot, and so we can color by state
#2015
transposed_scoring_data_2015$IsolateState <- c("ND",	"ND",	"ND",	"ND",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"SD",	"SD",	"SD",	"SD",	"NE",	"NE",	"NE",	"NE",	"OH",	"MS",	"TX",	"FL",	"FL")
#1990
transposed_scoring_data_1990$IsolateState <- c("MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"WI",	"WI",	"SD",	"SD",	"SD",	"PA",	"KS",	"AR",	"GA",	"LA",	"TX",	"TX",	"TX",	"TX",	"TX")

#Convert scoring data back to numeric
transposed_scoring_data_2015[,1:40] <- sapply(transposed_scoring_data_2015[,1:40], as.numeric)
transposed_scoring_data_1990[,1:40] <- sapply(transposed_scoring_data_1990[,1:40], as.numeric)

#Combine 1990 and 2015 for another PCA
transposed_scoring_data_1990_2015 <- rbind(transposed_scoring_data_2015, transposed_scoring_data_1990)

#Run PCA
#Use unscaled because variables all have the same units (Scoring data)
scorePCAmethods_2015 <- pca(transposed_scoring_data_2015[,1:40], nPcs = 2, method = "svd")
scorePCAmethods_1990 <- pca(transposed_scoring_data_1990[,1:40], nPcs = 2, method = "svd")
scorePCAmethods_1990_2015 <- pca(transposed_scoring_data_1990_2015[,1:40], nPcs = 2, method = "svd")
#str(scorePCAmethods)
#scorePCAmethods@R2

#Here's a quick way to make plots to examine side by side, not the best for saving
#slplot(scorePCAmethods, scoresLoadings = c(T,T), scol = transposed_scoring_data$StateColor, scex = 0.6, lcex = 0.6)

#Use ggplot
#First, merge by rows
df_2015 <- merge(scores(scorePCAmethods_2015), transposed_scoring_data_2015, by="row.names")
df_1990 <- merge(scores(scorePCAmethods_1990), transposed_scoring_data_1990, by="row.names")
df_1990_2015 <- merge(scores(scorePCAmethods_1990_2015), transposed_scoring_data_1990_2015, by="row.names")

library(ggplot2)

#Add rownames back
rownames(df_2015) <- c(df_2015$Row.names)
rownames(df_1990) <- c(df_1990$Row.names)
rownames(df_1990_2015) <- c(df_1990_2015$Row.names)

df_1990_2015$IsolateYear <- as.factor(c("2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015",  "2015", "2015", "2015", "2015", "2015", "2015", "2015", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990",  "1990", "1990", "1990", "1990", "1990", "1990", "1990"))

# Save the PCA plots to a file
pdf("PCA_plot_2015.pdf")

#Color by state
ggplot(df_2015, aes(x=df_2015$PC1, y=df_2015$PC2, color = df_2015$IsolateState)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("MN" = "firebrick2", "SD" = "gold", "TX" = "dodgerblue", "ND" = "black", "NE" = "darkorchid", "OH" = "darkgreen", "MS" = "darkorange", "FL" = "brown", "WI" = "darkgrey", "PA" = "deeppink", "KS" = "darkolivegreen4", "AR" = "black", "GA" = "tomato3", "LA" = "slateblue2")) +
  theme_classic(base_size = 16) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab(paste("PC1", scorePCAmethods_2015@R2[1] * 100, "% of variance")) +
  ylab(paste("PC2", scorePCAmethods_2015@R2[2] * 100, "% of variance"))

dev.off()

# Save the PCA plots to a file
pdf("PCA_plot_1990.pdf")

#Color by state
ggplot(df_1990, aes(x=df_1990$PC1, y=df_1990$PC2, color = df_1990$IsolateState)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("MN" = "firebrick2", "SD" = "gold", "TX" = "dodgerblue", "ND" = "black", "NE" = "darkorchid", "OH" = "darkgreen", "MS" = "darkorange", "FL" = "brown", "WI" = "darkgrey", "PA" = "deeppink", "KS" = "darkolivegreen4", "AR" = "black", "GA" = "tomato3", "LA" = "slateblue2")) +
  theme_classic(base_size = 16) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab(paste("PC1", scorePCAmethods_1990@R2[1] * 100, "% of variance")) +
  ylab(paste("PC2", scorePCAmethods_1990@R2[2] * 100, "% of variance"))

dev.off()

# Save the PCA plots to a file
pdf("PCA_plot_1990_2015.pdf")

#Color by state
ggplot(df_1990_2015, aes(x=df_1990_2015$PC1, y=df_1990_2015$PC2, color = df_1990_2015$IsolateYear)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("2015" = "firebrick2", "1990" = "black")) +
  theme_classic(base_size = 16) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab(paste("PC1", scorePCAmethods_1990_2015@R2[1] * 100, "% of variance")) +
  ylab(paste("PC2", scorePCAmethods_1990_2015@R2[2] * 100, "% of variance"))

dev.off()

# Save the loadings plots to a file

df_loading_2015 <- as.data.frame(loadings(scorePCAmethods_2015))
df_loading_1990 <- as.data.frame(loadings(scorePCAmethods_1990))
df_loading_1990_2015 <- as.data.frame(loadings(scorePCAmethods_1990_2015))

library(ggrepel)

pdf("PCA_loadings_plot_2015.pdf")

ggplot(df_loading_2015) +
  geom_point(aes(PC1, PC2), color = "grey", size = 3) +
  geom_text_repel(aes(PC1, PC2, label = rownames(df_loading_2015))) +
  theme_classic(base_size = 16)

dev.off()

pdf("PCA_loadings_plot_1990.pdf")

ggplot(df_loading_1990) +
  geom_point(aes(PC1, PC2), color = "grey", size = 3) +
  geom_text_repel(aes(PC1, PC2, label = rownames(df_loading_1990))) +
  theme_classic(base_size = 16)

dev.off()

pdf("PCA_loadings_plot_1990_2015.pdf")

ggplot(df_loading_1990_2015) +
  geom_point(aes(PC1, PC2), color = "grey", size = 3) +
  geom_text_repel(aes(PC1, PC2, label = rownames(df_loading_1990_2015))) +
  theme_classic(base_size = 16)

dev.off()

#PC1 explains most of the variation in the scoring data/phenotypes, PC2 the next highest amount, and so on
#Loadings are basically the influence scores of the differential lines that explain the variation along a direction

#Different plotting option
#library(ggrepel)
# ggplot(df) +
#   geom_point(aes(PC1, PC2), color = "grey", size = 3) +
#   geom_label_repel(aes(PC1, PC2, fill = IsolateState, label = rownames(df)),
#                    fontface = 'bold', color = 'white',
#                    box.padding = 0.35, point.padding = 0.5,
#                    segment.color = 'grey50') +
#   theme_classic(base_size = 16) +
#   scale_fill_manual(guide = guide_legend(override.aes = aes(color = NA)), values = setNames(c("darkorange3", "tomato3", "darkolivegreen4", "slateblue2", "firebrick2", "deeppink", "gold", "dodgerblue", "darkgrey"), levels(df$IsolateState))) +
#   xlab(paste("PC1", scorePCAmethods@R2[1] * 100, "% of variance")) +
#   ylab(paste("PC2", scorePCAmethods@R2[2] * 100, "% of variance")) +
#   labs(fill="Isolate State")