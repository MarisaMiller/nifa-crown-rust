#This short script is just for running pca on MSI to speed up computation

#This is for data analysis with adegenet using the information at http://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html, and also in the adegenet vignette.

#Only need vcfR for reading in original vcf files and doing some filtering
library("vcfR")
library("adegenet")
library("parallel")

load("gl_1990_rmNA.Rdata")
load("gl_2015_rmNA.Rdata")

pca_1990 <- glPca(gl_1990_rmNA, parallel = TRUE, n.cores = 15, nf=20)
save(pca_1990, file="pca_1990.Rdata")

pca_2015 <- glPca(gl_2015_rmNA, parallel = TRUE, n.cores = 15, nf=20)
save(pca_2015, file="pca_2015.Rdata")
