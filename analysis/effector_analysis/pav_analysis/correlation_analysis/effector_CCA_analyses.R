#############################################################################################################################
#############################################################################################################################
############################################### PHENOTYPE-GENOTYPE ASSOCIATIONS USING CCA ###################################
#############################################################################################################################
#############################################################################################################################
# Hoa Nguyen-Phuc and Marisa Miller

#############################################################################################################################
########################################################## PART 1 ###########################################################
################################################### PREPARING DATA FILES ####################################################
#############################################################################################################################

# Database management packages
library(gdata) # to read XLS files

# Genetics packages
library(ade4)      
library(adegenet)
library(ecodist)   # "distance"
library(ape)       # "as.dist"
library(phangorn)
library(phytools)
library(vcfR)
library(vegan) #For CCA

# Plotting
library(ggplot2)
library(ggvegan)

########################################################### DATA FILES ######################################################

### PHENOTYPES --------------------------------------------------------------------------------------------------------------
# First read in phenotyping data
# The column names correspond to the isolate names and the rows correspond to the different differential lines (Pc14, etc)
pheno1990 = as.matrix(read.csv("1990_phenotypes.csv", header=TRUE, check.names=FALSE, row.names = 1))
pheno2015 = as.matrix(read.csv("2015_phenotypes.csv", header=TRUE, check.names=FALSE, row.names = 1))

attributes(str(pheno1990))
attributes(str(pheno2015))

# Combine matrices for both years
phe = cbind(pheno1990,pheno2015)
dim(phe)

### PRESENCE/ABSENCE VARIATION -----------------------------------------------------------------------------------------------
# Read in presence/absence variation (coverage data) for haustorial effectors
# Columns correspond to the isolate names. Rows correspond to the gene names of haustorial effectors
# Values correspond to the proportion of the gene that had mapping coverage of less than 5: 0 (complete presence) to 1 (fully deleted)
pav = as.matrix(read.csv("effector_pav_data.csv", header=TRUE, check.names=FALSE, row.names = 1))
# Convert to a 0, 1 matrix and flip the values: 0 for any deletion and 1 for full presence
pav = ifelse(pav == 0, 1, 0)

dim(pav)

############################################### MATRIX & STRING MANIPULATIONS ###############################################

# Sort the data sets by column names
phe = phe[,order(colnames(phe))]
pav = pav[,order(colnames(pav))]

# Sort the variation set by row names
pav = pav[order(rownames(pav)),]

# Need to rotate the matrices 
# with Rows representing Samples (Isolates) and Columns representing Phenotypes or Genotypes
phe = t(phe)
pav = t(pav)

# Attributes - to make sure our three data sets have same dimensions and same col and row orders
dim(phe)
dim(pav)

colnames(phe)
colnames(pav)

row.names(phe)
row.names(pav)

# Rename col names in "pav"
gene.name.short = substr(colnames(pav),11,16)
gene.name.short
colnames(pav) = gene.name.short

# Assign population/cluster attributes
pop.year = substr(rownames(phe),1,2)
pop.state = substr(rownames(phe),3,4)

# Dealing with non-polymorphic --------------------------------------------------------------------------------------------

# "phe" data set
colSums(phe)
range(colSums(phe))
length(which((colSums(phe)==0)=="TRUE")) 
length(which((colSums(phe)==1)=="TRUE"))

# "pav" data set
# Remove non-polymorphic columns
colSums(pav)
length(which((colSums(pav)==0)=="TRUE")) 
length(which((colSums(pav)==nrow(pav))=="TRUE")) # 90 - these are non-polymorphic sites, need to remove them
pav = pav[,colSums(pav) != nrow(pav)]

setdiff(row.names(pav), row.names(phe))

#############################################################################################################################
########################################################## PART 2 ###########################################################
##################################################### DISTANCE MATRICES #####################################################
#############################################################################################################################

# Distance Class -----------------------------------------------------------------------------------------------------------
pav.dist = distance(pav, method="jaccard")
phe.dist = distance(phe, method="euclidean")

length(pav.dist)
length(phe.dist)

attributes(pav.dist)
attributes(phe.dist)

# Matrix Class -------------------------------------------------------------------------------------------------------------
pav.matrix = as.matrix(pav.dist)
phe.matrix = as.matrix(phe.dist)

dim(pav.matrix)
dim(phe.matrix)

head(pav.matrix)
head(phe.matrix)

# Optional Visualization -------------------------------------------------------------------------------------------------------------

# library(qgraph)
# 
# jpeg('phenotypes.jpg', width=3000, height=3000, unit='px')
# qgraph(1/phe.dist, layout='spring', vsize=3)
# dev.off()
# 
# jpeg('pre.abs.variation.jpg', width=3000, height=3000, unit='px')
# qgraph(1/pav.dist, layout='spring', vsize=3)
# dev.off()

#############################################################################################################################
########################################################## PART 3 ###########################################################
####################################################### ORDINATION ##########################################################
#############################################################################################################################

############################################## UNCONSTRAINSED ORDINATION ####################################################

# This part is for exploratory analysis using PCA
# dudi.pav.pca
# method: duality diagrm 
# nf is number of axes (principal components) retained

# Genoytopes ----------------------------------------------------------------------------------------------------------------

# # Full dataset of 2 vs 3 principal components
# pav.pca.dud = dudi.pca(pav, cent=FALSE, scannf=FALSE, nf=3)     
# pav.pca.prc = prcomp(pav)
# 
# # Describe all significant PCs
# attributes(str(pav.pca.dud))
# pav.pca.dud$eig
# pav.pca.dud$eig >1
# 
# # Principal Components
# pav.pca.dud$eig[1]/sum(pav.pca.dud$eig)
# pav.pca.dud$eig[2]/sum(pav.pca.dud$eig)
# pav.pca.dud$eig[3]/sum(pav.pca.dud$eig)
# sum(pav.pca.dud$eig[1:3])/sum(pav.pca.dud$eig)
# 
# # Plots With library "pav.pca3d"
# library(pca3d)
# library(rgl)
# # This has colored labels but can't change cex size
# # Create a matrix data, and change row names accordingly
# pav.pca3dmatrix = as.matrix(pav.pca.dud$l1)
# color3d = c("red","darkblue")
# # More customized plot (with defined col codes)
# pav.pca3Dframe2 = data.frame(pav.pca3dmatrix, factor(pop.year))
# attach(pav.pca3Dframe2); pca3d(pav.pca3dmatrix, col=color3d[factor(pop.year)], show.scale=TRUE,
#                                show.labels=rownames(pav), labels.col=color3d[factor(pop.year)])
# # wITHOUT LABELS                           
# attach(pav.pca3Dframe2); pca3d(pav.pca3dmatrix, col=color3d[factor(pop.year)], show.scale=TRUE)
# 
# # Phenoytopes ----------------------------------------------------------------------------------------------------------------
# phe.pca.dud = dudi.pca(phe, cent=FALSE, scannf=FALSE, nf=3)     
# 
# # Describe all significant PCs
# attributes(str(phe.pca.dud))
# phe.pca.dud$eig
# phe.pca.dud$eig >1
# 
# # Principal Components
# phe.pca.dud$eig[1]/sum(phe.pca.dud$eig)
# phe.pca.dud$eig[2]/sum(phe.pca.dud$eig)
# phe.pca.dud$eig[3]/sum(phe.pca.dud$eig)
# sum(phe.pca.dud$eig[1:3])/sum(phe.pca.dud$eig)
# 
# # Plots With library "phe.pca3d"
# # This has colored labels but can't change cex size
# # Create a matrix data, and change row names accordingly
# phe.pca3dmatrix = as.matrix(phe.pca.dud$l1)
# color3d = c("red","darkblue")
# # More customized plot (with defined col codes)
# phe.pca3Dframe2 = data.frame(phe.pca3dmatrix, factor(pop.year))
# attach(phe.pca3Dframe2); pca3d(phe.pca3dmatrix, col=color3d[factor(pop.year)], show.scale=TRUE,
#                                show.labels=rownames(phe), labels.col=color3d[factor(pop.year)])
# # wITHOUT LABELS                           
# attach(phe.pca3Dframe2); pca3d(phe.pca3dmatrix, col=color3d[factor(pop.year)], show.scale=TRUE)
# 
# # Eigenvalues ---------------------------------------------------------------------------------------------------------------
# par(mfrow=c(1,2))
# screeplot(pav.pca.dud)
# screeplot(phe.pca.dud)

#################################################### CONSTRAINSED ORDINATION ################################################

cca = cca(X=phe,Y=pav) # X is Response Variable, Y is Predictor
rda = rda(X=phe,Y=pav)

cca
rda

attributes(cca)
attributes(str(cca))

cca$inertia
cca$CCA

attributes(str(cca$CCA))
cca$CCA$eig
cca$CCA$u

dim(phe) # 60 x 40 # (40-1 = 39 degree of freedom)
dim(pav) # 60 x 12 # (12-1 = 11 degree of freedom)

cca$CCA$eig[1]/sum(cca$CCA$eig) # 44.9%
cca$CCA$eig[2]/sum(cca$CCA$eig) # 15.6%

cca$CA$eig[1]/sum(cca$CA$eig) # 20.0%
cca$CA$eig[2]/sum(cca$CA$eig) # 18.6%

# Scaling 1: Species scores (here are our Differential Lines) scaled to relative Eigenvalues
# Sites (here is our Samples) are weighted averages of the Species (here are our Differential Lines)
# Scaling 2: Site scores (here are our Samples) scaled to relative Eigenvalues
# Species (here are our Differential Lines) are weighted averages of the Sites

# an object found near the Centroid is likely to be found frequently

#Set theme to be classic
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=12),
              axis.title = element_text(size=15),
              strip.text.x = element_text(size = 13)
            )
)

cca_plot <- autoplot(cca, layers = c("species", "biplot"), geom = c("text")) + theme(legend.position="none")

pdf("cca.pdf")
cca_plot
dev.off()

#############################################################################################################################
#############################################################################################################################
########################################################### END #############################################################
#############################################################################################################################
#############################################################################################################################

#############################################################################################################################
####################################################### APPENDIX 1 ##########################################################
#############################################################################################################################

# TESTING ANOVA 

# anova(cca, step=1000)
# 
# par(mfrow=c(1,1))
# plot(cca, scaling=1, main = "Biplot CCA - scaling 1", display=c("sp","cn"),
#      xlab="", xaxt="n", ylab="", yaxt="n", ylim=c(0,0.8))
# #mtext("Geographic Distance (km)", side=1, line=3, cex=2)
# #mtext("Genetic Dissimilarity",    side=2, line=3, cex=2)
# axis(1, cex.axis=2)
# axis(2, cex.axis=2)
# 
# install.packages("devtools")
# devtools::install_github("gavinsimpson/ggvegan")
# 
# library(ggvegan)
# library(ggfortify)
# 
# autoplot(cca)
# autoplot(pav.pca.prc)
# 
# test = fortify(cca)
# test
# ggplot(test)
# 
# class(test)
# 
# autoplot(cca, colour = 'Species',
#          loadings = TRUE, loadings.colour = 'blue',
#          loadings.label = TRUE, loadings.label.size = 3)
# 
# autoplot(pav.pca.dud)

#############################################################################################################################
######################################################## APPENDIX 2 #########################################################
#############################################################################################################################

# TESTING MANTEL 

####################################################### ALL VECTORS #########################################################

# mantel = mantel.randtest(pav.dist, phe.dist, nrepet=99999)
# mantel
# 
# par(mfrow=c(1,1))
# plot(mantel, main="Mantel Test Phenotypes vs PAV", nclass=900)
# 
# corr = cor(phe.dist,pav.dist)
# corr
# which(corr == max(corr, na.rm = TRUE), arr.ind = TRUE)
# 
# ############################################### SINGLE-VECTOR DISTANCE MATRIX ##############################################
# 
# pav1 = distance(pav[1,]) # 102 rows for 102 genes
# phe1 = distance(phe[1,])     #  40 rows for  40 diffential lines   
# 
# length(pav1)
# length(phe1)
# 
# man1 = mantel.randtest(pav1, phe1, nrepet=99999)
# man1
# 
# lapply(1:length(pair), function(i) head(pair[[i]]))
# 
# pav.vector = lapply(1:nrow(pav), function(i) distance(pav[i,]))
# phe.vector = lapply(1:nrow(phe), function(i) distance(phe[i,]))
# 
# man.vector.phe = lapply(1:nrow(phe), function(i) mantel.randtest(pav.vector[[2]], phe.vector[[i]], nrepet=99999))
# man.vector.phe

#############################################################################################################################
#############################################################################################################################
########################################################### END #############################################################
#############################################################################################################################
#############################################################################################################################
























