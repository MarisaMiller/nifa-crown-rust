#This is for data analysis with adegenet using the information at http://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html, and also in the adegenet vignette.

#Change working directory to appripriate location
#setwd("~/Documents/Work/Log Book Data/2016/2016_crownrust_NIFA_project/sequencing_data_analysis/poppr_adegenet_dapc_analysis")

#Only need vcfR for reading in original vcf files and doing some filtering
#library("vcfR")
library("adegenet")
library("parallel")
#library("gridExtra") #2.3
#library("gridGraphics") #0.2
library("ape") #5.0
#library("ggtree") #1.8.2
#library("matrixStats")

load("gl_1990_2015_rmNA.Rdata")

#keep 5 principal components
#pca_1990_2015 <- glPca(gl_1990_2015_rmNA, parallel = TRUE, n.cores = 20, nf=20)

#save(pca_1990, file="pca_1990.Rdata")
#save(pca_2015, file="pca_2015.Rdata")
#save(pca_1990_2015, file="pca_1990_2015.Rdata")

dist <- dist(as.matrix(gl_1990_2015_rmNA))
tree_1990_2015 <- nj(dist)

tree_1990_2015_boots <- boot.phylo(tree_1990_2015, as.matrix(gl_1990_2015_rmNA), function(e) root(nj(dist(e)), 1), B = 1000, mc.cores = 16)
save(tree_1990_2015_boots, file="tree_1990_2015_boots.Rdata")


