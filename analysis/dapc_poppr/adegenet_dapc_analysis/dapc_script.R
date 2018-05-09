#This is for data analysis with adegenet using the information at http://grunwaldlab.github.io/Population_Genetics_in_R/analysis_of_genome.html and also in the adegenet tutorials here (https://github.com/thibautjombart/adegenet/wiki/Tutorials).

#This is a very interactive analysis, read script comments carefully.

#Change working directory to appropriate location
setwd("")

#Only need vcfR for reading in original vcf files and doing some filtering
library("vcfR") #1.8.0
library("adegenet") #2.1.1
library("gridExtra") #2.3
library("gridGraphics") #0.2
library("ape") #5.0
library("ggtree") #1.8.2
library("matrixStats")

#Read in VCF files for 1990 and 2015 isolates, just do once and then convert to genlight (see below)
#The combined year vcf was created with vcfcombine in vcflib
vcf_1990 <- read.vcfR("../../snp_calling/freebayes/1990_isolates.filter.vcf")
vcf_2015 <- read.vcfR("../../snp_calling/freebayes/2015_isolates.filter.vcf")
vcf_1990_2015 <- read.vcfR("../../snp_calling/freebayes/1990_2015_isolates.filter.vcf")

#Convert the large and unwieldy VCF files into smaller genlight objects that can be used with adegenet
#Non-biallelic variants will be discarded since genlight objects don't support non-biallelic data
gl_1990 <- vcfR2genlight(vcf_1990)
gl_2015 <- vcfR2genlight(vcf_2015)
gl_1990_2015 <- vcfR2genlight(vcf_1990_2015)

save(gl_1990, file="gl_1990.Rdata")
save(gl_2015, file="gl_2015.Rdata")
save(gl_1990_2015, file="gl_1990_2015.Rdata")

#Then reload objects for subsequent analyses if need be
# load("gl_1990.Rdata")
# load("gl_2015.Rdata")
# load("gl_1990_2015.Rdata")

#Set population information to samples
#This will be a case-by-case basis for your samples
pop(gl_1990) <- as.factor(c("GA", "TX", "TX", "MN-B", "MN-B", "MN", "MN", "MN-B", "KS", "MN", "WI", "MN-B", "MN-B", "TX", "SD", "MN-B",  "WI", "MN-B", "MN-B", "SD", "MN", "MN-B", "LA", "MN", "PA", "TX", "TX", "SD", "AR", "MN-B"))

pop(gl_2015) <- as.factor(c("ND", "MN", "NE", "ND", "MN", "NE", "FL", "MN-B", "SD", "MN-B", "SD", "MN-B", "SD", "NE", "ND", "MN", "OH", "TX", "SD", "FL", "NE", "MN-B", "MN", "MN-B", "MN", "12SD80", "MN-B", "ND", "12NC29", "MS", "MN", "MN-B"))

pop(gl_1990_2015) <- as.factor(c("ND", "MN", "NE", "ND", "MN", "NE", "FL", "MN-B", "SD", "MN-B", "SD", "MN-B", "SD", "NE", "ND", "MN", "OH", "TX", "SD", "FL", "NE", "MN-B", "MN", "MN-B", "MN", "12SD80", "MN-B", "ND", "12NC29", "MS", "MN", "MN-B", "GA", "TX", "TX", "MN-B", "MN-B", "MN", "MN", "MN-B", "KS", "MN", "WI", "MN-B", "MN-B", "TX", "SD", "MN-B",  "WI", "MN-B", "MN-B", "SD", "MN", "MN-B", "LA", "MN", "PA", "TX", "TX", "SD", "AR", "MN-B"))

#Remove reference isolates (12NC29 and 12SD80) for this analysis and drop alleles no longer polymorphic in the data
gl_2015 <- gl_2015[!(indNames(gl_2015) %in% c("12SD80", "12NC29"))]
gl_1990_2015 <- gl_1990_2015[!(indNames(gl_1990_2015) %in% c("12SD80", "12NC29"))]

#Set ploidy
gl_list <- list(gl_1990, gl_2015, gl_1990_2015)
lapply(gl_list, ploidy, 2)

#First remove any NAs
toRemove1990 <- is.na(glMean(gl_1990, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove1990) # position of entirely non-typed loci
gl_1990_rmNA <- gl_1990[, !toRemove1990]

toRemove2015 <- is.na(glMean(gl_2015, alleleAsUnit = FALSE))
which(toRemove2015)
gl_2015_rmNA <- gl_2015[, !toRemove2015]

toRemove1990_2015 <- is.na(glMean(gl_1990_2015, alleleAsUnit = FALSE))
which(toRemove1990_2015)
gl_1990_2015_rmNA <- gl_1990_2015[, !toRemove1990_2015]

save(gl_1990_rmNA, file="gl_1990_rmNA.Rdata")
save(gl_2015_rmNA, file="gl_2015_rmNA.Rdata")
save(gl_1990_2015_rmNA, file="gl_1990_2015_rmNA.Rdata")

#Then reload objects for subsequent analyses if need be
load("gl_1990_rmNA.Rdata") #1,379,292 biallelic SNPs used for analysis
load("gl_2015_rmNA.Rdata") #1,186,114 biallelic SNPs
load("gl_1990_2015_rmNA.Rdata") #2,007,463 biallelic SNPs

#Let's add some strata to the gl object, then we can be more flexible in how we define our populations for DAPC.
#BT means BT in more than 3 counties (according to here: https://www.eddmaps.org/distribution/usstate.cfm?sub=3070), accessed Feb 9 2018
#Citation: EDDMapS. 2018. Early Detection & Distribution Mapping System.
#For 1990 isolates, we aren't sure if the BT distribution was similar to now, but we can assume it is so (https://link.springer.com/article/10.1007/s10530-007-9091-3).
strata(gl_1990_rmNA) <- data.frame(BT=c("NO BT", "NO BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "BT", "NO BT", "BT"), BT2=c("NO BT", "NO BT", "NO BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "MN BT", "MN BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "BT", "NO BT", "MN BT"))

strata(gl_2015_rmNA) <- data.frame(BT=c("BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT"), BT2=c("BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "NO BT", "BT", "MN BT"))

strata(gl_1990_2015_rmNA) <- data.frame(BT=c("BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "BT", "NO BT", "BT"), BT2=c("BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "NO BT", "BT", "MN BT", "NO BT", "NO BT", "NO BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "MN BT", "MN BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "BT", "NO BT", "MN BT"), year=c(rep("2015", 30), rep("1990", 30)))

#Perform PCA in adegenet to speed up computation for finding clusters for DAPC
#From the adegent manual "The trade-off between power of discrimination and over-fitting can be measured by the a-score, which is simply the difference between the proportion of successful reassignment of the analysis (observed discrimination) and values obtained using random groups (random discrimination)."
#"The a-score can serve as a criterion for choosing the optimal number of PCs in the PCA step of DAPC, i.e. the number of PC maximizing the a-score."
#I used the optim.a.score function in adegenet to compute the a-score. I used the default setting (employing the "smart" procedure to predict the approximate optimal number of PCs." 
#I explored the a-score plot (see below) and identified the ideal number of PCA components to keep for each analysis.
pca_1990 <- glPca(gl_1990_rmNA, useC = TRUE, nf=20)
pca_2015 <- glPca(gl_2015_rmNA, useC = TRUE, nf=20)
pca_1990_2015 <- glPca(gl_1990_2015_rmNA, useC = TRUE, nf=20)

#The PCA can take a while, so run on a cluster with more compute power if possible
#Here's an example
#library("parallel")
#pca_2015 <- glPca(gl_2015_rmNA, parallel = TRUE, n.cores = 15, nf=20)

save(pca_1990, file="pca_1990.Rdata")
save(pca_2015, file="pca_2015.Rdata")
save(pca_1990_2015, file="pca_1990_2015.Rdata")

#Load for subsequent analysis if need be
load("pca_1990.Rdata")
load("pca_2015.Rdata")
load("pca_1990_2015.Rdata")

#Perform DAPC in adegenet
#Set pop by BT, MN-BT (where we know for sure there is BT), and NON-BT, rather than by just states, since these are our biologically informative "groups"
setPop(gl_1990_rmNA) <- ~BT2
setPop(gl_2015_rmNA) <- ~BT2
setPop(gl_1990_2015_rmNA) <- ~BT2

#Visualize gl SNP data results for combined years as neighbor joining tree
dist <- dist(as.matrix(gl_1990_2015_rmNA))
tree_1990_2015 <- nj(dist)

#Is NJ a good fit for our data? To check, follow the suggestions of Jombart in this tutorial: http://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf
#This is easily investigated using simple biplots and correlation indices.
#The function cophenetic is used to compute distances between the tips of the tree. 
x <- as.vector(dist)
y <- as.vector(as.dist(cophenetic(tree_1990_2015)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", main="Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2 #0.9887178, so it is a good fit

#Bootstrap with 1000 replicates
#This can be computationally expensive, so save data afterwards to reload if need be
#Load the parallel library if using the mc.cores option
#library(parallel)
tree_1990_2015_boots <- boot.phylo(tree_1990_2015, as.matrix(gl_1990_2015_rmNA), function(e) root(nj(dist(e)), 1), B = 1000, mc.cores = 16)
save(tree_1990_2015_boots, file="tree_1990_2015_boots.Rdata")
load("tree_1990_2015_boots.Rdata")

#Convert to percent
tree_1990_2015_boots <- (tree_1990_2015_boots/1000)*100
#Filter out nodes less than 75% boostrap support
tree_1990_2015_boots <- unlist(lapply(tree_1990_2015_boots, function(x) {if(x>75){x<-NA} else (x)}))

#Add bootstraps to tree
tree_1990_2015$node.label <- tree_1990_2015_boots

#Add bootstraps to treedata
#tree2treedata <- treeio::as.treedata(tree_1990_2015, tree_1990_2015_boots)

popColorInfo <- data.frame(isolate = tree_1990_2015$tip.label, population = gl_1990_2015_rmNA@pop, color = NA)
popColorInfo$color[popColorInfo$population == "NO BT"] <- "#E69F00"
popColorInfo$color[popColorInfo$population == "BT"] <- "#D55E00"
popColorInfo$color[popColorInfo$population == "MN BT"] <- "#56B4E9"

groupInfo <- split(tree_1990_2015$tip.label, gsub("([1590]{2}).+", "\\1", tree_1990_2015$tip.label))
tree_1990_2015 <- groupOTU(tree_1990_2015, groupInfo)

tree_plot <- ggtree(tree_1990_2015, aes(color = group), layout='circular') + geom_tiplab(size=3, aes(angle=angle), color = "black", offset = 20) + geom_point2(aes(label=label, subset=!isTip & as.numeric(label) < 80), color = "green", size=1)

tree_plot <- tree_plot %<+% popColorInfo + geom_tippoint(aes(color = color), size=2) + scale_color_manual(values = c("#56B4E9", "#D55E00", "#E69F00", "red", "black")) + geom_treescale(x=1000, y=0, width = 100) #The scale represents nucleotide differences

#Save plot
pdf("nj_tree_1990_2015.pdf", height = 7, width = 7)
tree_plot
dev.off()

#Next, run dapc with all PCs retained to examine a-score at each level. 
#Keep all discriminant functions, then examine plot to determine how many PCs to use
dapc_1990 <- dapc(gl_1990_rmNA, gl_1990_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_1990)
temp_1990 <- optim.a.score(dapc_1990) #Use 10 PCs since the spline interpolation suggested that as the ideal amount.
#Now, re-run Discriminant Analysis on the retained principal components, with the appropriate amount suggested by optim.a.score 
dapc_1990 <- dapc(gl_1990_rmNA, gl_1990_rmNA@pop, n.pca = 10, n.da = 3, glPca = pca_1990)

dapc_2015 <- dapc(gl_2015_rmNA, gl_2015_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_2015)
temp_2015 <- optim.a.score(dapc_2015) #Use 9 PCs
dapc_2015 <- dapc(gl_2015_rmNA, gl_2015_rmNA@pop, n.pca = 9, n.da = 3, glPca = pca_2015)

dapc_1990_2015 <- dapc(gl_1990_2015_rmNA, gl_1990_2015_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_1990_2015)
temp_1990_2015 <- optim.a.score(dapc_1990_2015) #Use 15 PCs
dapc_1990_2015 <- dapc(gl_1990_2015_rmNA, gl_1990_2015_rmNA@pop, n.pca = 15, n.da = 3, glPca = pca_1990_2015)

#Make plots of results
##1990
#The center of each group is indicated with crosses.

#Make a function to save plots
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

myCol = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")
par(bty = "n") #Remove box around plots
scatter(dapc_1990, lwd = 2, lty = 2, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = myCol, xlim = c(-4,4))
points(dapc_1990$grp.coord[,1], dapc_1990$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_1990$grp.coord[,1], dapc_1990$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=myCol)
scatter_1990 <- grab_grob()

#Reformat so that individuals are sorted by their expected membership so that compoplot is easier to interpret
#These plots show membership probabilities of each individual for the different biological groups based on the retained discriminant functions
#Membership probabilities also provide indications of how clear-cut genetic clusters are.
#Loose clusters will result in fairly flat distributions of membership probabilities of individuals across clusters, pointing to possible admixture.
compoplot_1990_df <- data.frame(grp=dapc_1990$grp, dapc_1990$posterior)
compoplot_1990_df <- compoplot_1990_df[order(compoplot_1990_df$grp),]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_1990_df[,2:4]), cex.axis = 0.75, cex.names = 0.5, col.pal = myCol, show.lab = TRUE, legend = FALSE)
compo_1990 <- grab_grob()

#Calculate average and standard deviation of admixture per group
colMeans(compoplot_1990_df[compoplot_1990_df$grp == "NO BT", 2:ncol(compoplot_1990_df)])
colMeans(compoplot_1990_df[compoplot_1990_df$grp == "BT", 2:ncol(compoplot_1990_df)])
colMeans(compoplot_1990_df[compoplot_1990_df$grp == "MN BT", 2:ncol(compoplot_1990_df)])

colSds(as.matrix(compoplot_1990_df[compoplot_1990_df$grp == "NO BT", 2:ncol(compoplot_1990_df)]))
colSds(as.matrix(compoplot_1990_df[compoplot_1990_df$grp == "BT", 2:ncol(compoplot_1990_df)]))
colSds(as.matrix(compoplot_1990_df[compoplot_1990_df$grp == "MN BT", 2:ncol(compoplot_1990_df)]))
#NO BT Isolates
# NO.BT     MN.BT        BT 
# Average
# 0.7177748 0.1703482 0.1118770
# SD
# 0.3412310 0.2605177 0.2959584

#BT Isolates
# NO.BT     MN.BT        BT 
# 0.1270517 0.0127306 0.8602177 
# SD
# 0.26640001 0.02445583 0.29068692

#MN BT Isolates
# NO.BT     MN.BT        BT 
# 0.23638877 0.69626019 0.06735104 
# SD
# 0.3430090 0.4152446 0.2010942

#Can also use assignplots to look at group memberships after DAPC. While not great for publication, can be useful to examine.
#Large values indicate clear-cut clusters, while low values suggest admixed (or even misclassified) groups.
#Heat colors represent membership probabilities (red=1, white=0); blue crosses represent the prior cluster provided to DAPC.
assignplot(dapc_1990)

##2015
#I have to use a different color palette since the order of grp levels is different and not worth re-ordering here (aka, I spent too long fiddling with it to no avail)
scatter(dapc_2015, lwd = 2, lty = 2, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = c("#D55E00", "#E69F00", "#56B4E9"))
legend(3, 2.5, c("NO BT", "MN BT", "BT"), bty = "n", pch = 20, col = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00"))
points(dapc_2015$grp.coord[,1], dapc_2015$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_2015$grp.coord[,1], dapc_2015$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=c("#D55E00", "#E69F00", "#56B4E9"))
scatter_2015 <- grab_grob()

compoplot_2015_df <- data.frame(grp=dapc_2015$grp, dapc_2015$posterior)
#Re-order factor levels and columns to be same order as 1990
compoplot_2015_df$grp = factor(compoplot_2015_df$grp,levels(compoplot_2015_df$grp)[c(2,3,1)])
compoplot_2015_df <- compoplot_2015_df[order(compoplot_2015_df$grp),]
compoplot_2015_df <- compoplot_2015_df[,c(1,3,4,2)]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_2015_df[,2:4]), col.pal = myCol, cex.axis = 0.75, cex.names = 0.5, show.lab = TRUE, legend = FALSE)
compo_2015 <- grab_grob()

colMeans(compoplot_2015_df[compoplot_2015_df$grp == "NO BT", 2:ncol(compoplot_2015_df)])
colMeans(compoplot_2015_df[compoplot_2015_df$grp == "BT", 2:ncol(compoplot_2015_df)])
colMeans(compoplot_2015_df[compoplot_2015_df$grp == "MN BT", 2:ncol(compoplot_2015_df)])

colSds(as.matrix(compoplot_2015_df[compoplot_2015_df$grp == "NO BT", 2:ncol(compoplot_2015_df)]))
colSds(as.matrix(compoplot_2015_df[compoplot_2015_df$grp == "BT", 2:ncol(compoplot_2015_df)]))
colSds(as.matrix(compoplot_2015_df[compoplot_2015_df$grp == "MN BT", 2:ncol(compoplot_2015_df)]))
#NO BT Isolates
# NO.BT     MN.BT        BT 
# 0.53234418 0.01942991 0.44822591
# SD
# 0.44639901 0.01864153 0.42785531

#BT Isolates
# NO.BT     MN.BT        BT 
# 0.09324644 0.10149686 0.80525670 
# SD
# 0.1988935 0.2178652 0.2744180

#MN BT Isolates
# NO.BT     MN.BT        BT 
# 0.02972409 0.64601995 0.32425597 
# SD
# 0.05413863 0.32950100 0.28389797

assignplot(dapc_2015)

#Save the plots
#The y-axis label for the compoplots has a bug with grid.arrange and is not displayed, so the label was added on in Adobe Illustrator
#Labels for the group membership of the individuals in the compoplot was also added on in Illustrator
pdf("1990_2015_dapc.pdf", height = 9, width = 7)
grid.arrange(arrangeGrob(scatter_1990, top = "1990", clip = "on"), arrangeGrob(scatter_2015, top = "2015", clip = "on"), arrangeGrob(compo_1990, clip = "on"), arrangeGrob(compo_2015, clip = "on"), ncol=2, nrow=2, vp = viewport(width=0.98, height=0.98))
dev.off()

#Plots for 1990_2015 combined data
scatter(dapc_1990_2015, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = c("#D55E00", "#E69F00", "#56B4E9"))
legend(-3, 2.5, c("NO BT", "MN BT", "BT"), bty = "n", pch = 20, col = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00"))
points(dapc_1990_2015$grp.coord[,1], dapc_1990_2015$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_1990_2015$grp.coord[,1], dapc_1990_2015$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=c("#D55E00", "#E69F00", "#56B4E9"))
scatter_1990_2015 <- grab_grob()

compoplot_1990_2015_df <- data.frame(grp=dapc_1990_2015$grp, dapc_1990_2015$posterior)
#Re-order factor levels and columns to be same order as 1990
compoplot_1990_2015_df$grp = factor(compoplot_1990_2015_df$grp,levels(compoplot_1990_2015_df$grp)[c(2,3,1)])
compoplot_1990_2015_df <- compoplot_1990_2015_df[order(compoplot_1990_2015_df$grp),]
compoplot_1990_2015_df <- compoplot_1990_2015_df[,c(1,3,4,2)]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_1990_2015_df[,2:4]), col.pal = myCol, cex.axis = 0.75, cex.names = 0.35, show.lab = TRUE, legend = FALSE)
compo_1990_2015 <- grab_grob()

pdf("1990_2015_combined_dapc_by_BT.pdf", height = 9, width = 5)
grid.arrange(arrangeGrob(scatter_1990_2015, top = "1990 and 2015 combined", clip = "on"), arrangeGrob(compo_1990_2015, clip = "on"), ncol=1, nrow=2, vp = viewport(width=0.95, height=0.95))
dev.off()

#Just for fun, we'll look at the combined dataset by year to see if years can be differentiated.
setPop(gl_1990_2015_rmNA) <- ~year #re-run the dapc command after re-assigning the pop
dapc_1990_2015 <- dapc(gl_1990_2015_rmNA, gl_1990_2015_rmNA@pop, n.pca = 15, n.da = 3, glPca = pca_1990_2015)
#Since there are only two groups for this data (1990 and 2015), the results are visualized by plotting the densities of individuals on a given discriminant function with different colors for different groups.
scatter(dapc_1990_2015, scree.da = FALSE, legend = FALSE, cex.axis = 0.75, col = c("2015" = "red", "1990" = "black"))
legend(-6, 1.75, c("2015", "1990"), bty = "n", pch = 20, col = c("2015" = "red", "1990" = "black"))
scatter_1990_2015 <- grab_grob()

par(mai=c(1.5,1,1,1))
compoplot_1990_2015_df <- data.frame(grp=dapc_1990_2015$grp, dapc_1990_2015$posterior)
compoplot_1990_2015_df <- compoplot_1990_2015_df[order(compoplot_1990_2015_df$grp),]
compoplot(as.matrix(compoplot_1990_2015_df[,2:3]), col.pal = c("2015" = "red", "1990" = "black"), cex.axis = 0.75, cex.names = 0.35, show.lab = TRUE, legend = FALSE)
compo_1990_2015 <- grab_grob()

pdf("1990_2015_combined_dapc_by_year.pdf", height = 9, width = 5)
grid.arrange(arrangeGrob(scatter_1990_2015, top = "1990 and 2015 combined", clip = "on"), arrangeGrob(compo_1990_2015, clip = "on"), ncol=1, nrow=2, vp = viewport(width=0.95, height=0.95))
dev.off()

#####Repeat analysis with VCF files that are the resulting of merging results from 3 variant callers for 1990 and 2015#####
vcf_merge_1990 <- read.vcfR("../../snp_calling/compare_gatk_freebayes_samtools/1990_gatk_freebayes_samtools_isolates.intersect.final.vcf")
vcf_merge_2015 <- read.vcfR("../../snp_calling/compare_gatk_freebayes_samtools/2015_gatk_freebayes_samtools_isolates.intersect.final.vcf")

#Convert to genlight objects 
gl_1990_merge <- vcfR2genlight(vcf_merge_1990)
gl_2015_merge <- vcfR2genlight(vcf_merge_2015)

#Set population information to samples
pop(gl_1990_merge) <- as.factor(c("AR", "GA", "KS", "LA", "MN", "MN-B", "MN", "MN", "MN-B", "MN", "MN", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "PA", "SD", "SD", "SD", "TX", "TX", "TX", "TX", "TX", "WI", "WI"))

pop(gl_2015_merge) <- as.factor(c("NC", "SD", "FL", "FL", "MN", "MN", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "MN-B", "MN", "MN", "MN", "MN", "MS", "ND", "ND", "ND", "ND", "NE", "NE", "NE", "NE", "OH", "SD", "SD", "SD", "SD", "TX"))

#Remove reference isolates (12NC29 and 12SD80) for this analysis and drop alleles no longer polymorphic in the data
gl_2015_merge <- gl_2015_merge[!(indNames(gl_2015_merge) %in% c("12SD80", "12NC29"))]

#Set ploidy
gl_merge_list <- list(gl_1990_merge, gl_2015_merge)
lapply(gl_merge_list, ploidy, 2)

#First remove any NAs
toRemove1990_merge <- is.na(glMean(gl_1990_merge, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove1990_merge) # position of entirely non-typed loci
gl_1990_merge_rmNA <- gl_1990_merge[, !toRemove1990_merge]

toRemove2015_merge <- is.na(glMean(gl_2015_merge, alleleAsUnit = FALSE))
which(toRemove2015_merge)
gl_2015_merge_rmNA <- gl_2015_merge[, !toRemove2015_merge]

save(gl_1990_merge_rmNA, file="gl_1990_merge_rmNA.Rdata")
save(gl_2015_merge_rmNA, file="gl_2015_merge_rmNA.Rdata")

#Then reload objects for subsequent analyses if need be
#load("gl_1990_merge_rmNA.Rdata") #2,008,276 biallelic SNPs used for analysis
#load("gl_2015_merge_rmNA.Rdata") #1,713,750 biallelic SNPs

#Let's add some strata to the gl object, then we can be more flexible in how we define our populations for DAPC.
strata(gl_1990_merge_rmNA) <- data.frame(BT=c("NO BT", "NO BT", "NO BT", "NO BT", "BT", "MN BT", "BT", "BT", "MN BT", "BT", "BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "BT", "BT", "BT", "BT", "NO BT", "NO BT", "NO BT", "NO BT", "NO BT", "BT", "BT"))

strata(gl_2015_merge_rmNA) <- data.frame(BT=c("NO BT", "NO BT", "BT", "BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "MN BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "BT", "NO BT"))

#Perform PCA in adegenet to speed up computation for finding clusters for DAPC
pca_merge_1990 <- glPca(gl_1990_merge_rmNA, useC = TRUE, nf=20)
pca_merge_2015 <- glPca(gl_2015_merge_rmNA, useC = TRUE, nf=20)

save(pca_merge_1990, file="pca_merge_1990.Rdata")
save(pca_merge_2015, file="pca_merge_2015.Rdata")

#Load for subsequent analysis if need be
#load("pca_merge_1990.Rdata")
#load("pca_merge_2015.Rdata")

#Perform DAPC in adegenet
#Set pop by BT, MN-BT (where we know for sure there is BT), and NON-BT, rather than by just states, since these are our biologically informative "groups"
setPop(gl_1990_merge_rmNA) <- ~BT
setPop(gl_2015_merge_rmNA) <- ~BT

#Next, run dapc with all PCs retained to examine a-score at each level. 
dapc_merge_1990 <- dapc(gl_1990_merge_rmNA, gl_1990_merge_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_merge_1990)
temp_merge_1990 <- optim.a.score(dapc_merge_1990) #Use 10 PCs 
dapc_merge_1990 <- dapc(gl_1990_merge_rmNA, gl_1990_merge_rmNA@pop, n.pca = 10, n.da = 3, glPca = pca_merge_1990)

dapc_merge_2015 <- dapc(gl_2015_merge_rmNA, gl_2015_merge_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_merge_2015)
temp_merge_2015 <- optim.a.score(dapc_merge_2015) #Use 10 PCs
dapc_merge_2015 <- dapc(gl_2015_merge_rmNA, gl_2015_merge_rmNA@pop, n.pca = 10, n.da = 3, glPca = pca_merge_2015)

#Make plots of results
##1990

#Make a function to save plots
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

myCol = c("NO BT" = "#E69F00", "BT" = "#D55E00", "MN BT" = "#56B4E9")
par(bty = "n") #Remove box around plots
scatter(dapc_merge_1990, lwd = 2, lty = 2, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = myCol, xlim = c(-4,4))
points(dapc_merge_1990$grp.coord[,1], dapc_merge_1990$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_merge_1990$grp.coord[,1], dapc_merge_1990$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=myCol)
scatter_merge_1990 <- grab_grob()

#Reformat so that individuals are sorted by their expected membership so that compoplot is easier to interpret
compoplot_merge_1990_df <- data.frame(grp=dapc_merge_1990$grp, dapc_merge_1990$posterior)
compoplot_merge_1990_df$grp = factor(compoplot_merge_1990_df$grp,levels(compoplot_merge_1990_df$grp)[c(1,3,2)])
compoplot_merge_1990_df <- compoplot_merge_1990_df[order(compoplot_merge_1990_df$grp),]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_merge_1990_df[,2:4]), cex.axis = 0.75, cex.names = 0.5, col.pal = myCol, show.lab = TRUE, legend = FALSE)
compo_merge_1990 <- grab_grob()

##2015
scatter(dapc_merge_2015, lwd = 2, lty = 2, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = c("NO BT" = "#E69F00", "BT" = "#D55E00", "MN BT" = "#56B4E9"))
legend(3, 2.5, c("NO BT", "MN BT", "BT"), bty = "n", pch = 20, col = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00"))
points(dapc_merge_2015$grp.coord[,1], dapc_merge_2015$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_merge_2015$grp.coord[,1], dapc_merge_2015$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=c("NO BT" = "#E69F00", "BT" = "#D55E00", "MN BT" = "#56B4E9"))
scatter_merge_2015 <- grab_grob()

compoplot_merge_2015_df <- data.frame(grp=dapc_merge_2015$grp, dapc_merge_2015$posterior)
#Re-order factor levels and columns to be same order as 1990
compoplot_merge_2015_df$grp = factor(compoplot_merge_2015_df$grp,levels(compoplot_merge_2015_df$grp)[c(1,3,2)])
compoplot_merge_2015_df <- compoplot_merge_2015_df[order(compoplot_merge_2015_df$grp),]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_merge_2015_df[,2:4]), col.pal = myCol, cex.axis = 0.75, cex.names = 0.5, show.lab = TRUE, legend = FALSE)
compo_merge_2015 <- grab_grob()

#Save the plots
pdf("1990_2015_merge_dapc.pdf", height = 9, width = 7)
grid.arrange(arrangeGrob(scatter_merge_1990, top = "1990", clip = "on"), arrangeGrob(scatter_merge_2015, top = "2015", clip = "on"), arrangeGrob(compo_merge_1990, clip = "on"), arrangeGrob(compo_merge_2015, clip = "on"), ncol=2, nrow=2, vp = viewport(width=0.98, height=0.98))
dev.off()

#####Analysis for just synonymous substitutions using FreeBayes only set#####
vcf_syn_1990 <- read.vcfR("../../snp_calling/variant_annotation/1990_isolates.filter.syn.vcf")
vcf_syn_2015 <- read.vcfR("../../snp_calling/variant_annotation/2015_isolates.filter.syn.vcf")

#Convert to genlight objects 
gl_1990_syn <- vcfR2genlight(vcf_syn_1990)
gl_2015_syn <- vcfR2genlight(vcf_syn_2015)

#Set population information to samples
pop(gl_1990_syn) <- as.factor(c("GA", "TX", "TX", "MN-B", "MN-B", "MN", "MN", "MN-B", "KS", "MN", "WI", "MN-B", "MN-B", "TX", "SD", "MN-B",  "WI", "MN-B", "MN-B", "SD", "MN", "MN-B", "LA", "MN", "PA", "TX", "TX", "SD", "AR", "MN-B"))

pop(gl_2015_syn) <- as.factor(c("ND", "MN", "NE", "ND", "MN", "NE", "FL", "MN-B", "SD", "MN-B", "SD", "MN-B", "SD", "NE", "ND", "MN", "OH", "TX", "SD", "FL", "NE", "MN-B", "MN", "MN-B", "MN", "12SD80", "MN-B", "ND", "12NC29", "MS", "MN", "MN-B"))

#Remove reference isolates (12NC29 and 12SD80) for this analysis and drop alleles no longer polymorphic in the data
gl_2015_syn <- gl_2015_syn[!(indNames(gl_2015_syn) %in% c("12SD80", "12NC29"))]

#Set ploidy
gl_syn_list <- list(gl_1990_syn, gl_2015_syn)
lapply(gl_syn_list, ploidy, 2)

#First remove any NAs
toRemove1990_syn <- is.na(glMean(gl_1990_syn, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove1990_syn) # position of entirely non-typed loci
gl_1990_syn_rmNA <- gl_1990_syn[, !toRemove1990_syn]

toRemove2015_syn <- is.na(glMean(gl_2015_syn, alleleAsUnit = FALSE))
which(toRemove2015_syn)
gl_2015_syn_rmNA <- gl_2015_syn[, !toRemove2015_syn]

save(gl_1990_syn_rmNA, file="gl_1990_syn_rmNA.Rdata")
save(gl_2015_syn_rmNA, file="gl_2015_syn_rmNA.Rdata")

#Then reload objects for subsequent analyses if need be
#load("gl_1990_syn_rmNA.Rdata") #74,976 biallelic SNPs used for analysis
#load("gl_2015_syn_rmNA.Rdata") #65,247 biallelic SNPs

#Let's add some strata to the gl object, then we can be more flexible in how we define our populations for DAPC.
strata(gl_1990_syn_rmNA) <- data.frame(BT=c("NO BT", "NO BT", "NO BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "MN BT", "MN BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "BT", "NO BT", "MN BT"))

strata(gl_2015_syn_rmNA) <- data.frame(BT=c("BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "NO BT", "BT", "MN BT"))

#Perform PCA in adegenet to speed up computation for finding clusters for DAPC
pca_syn_1990 <- glPca(gl_1990_syn_rmNA, useC = TRUE, nf=20)
pca_syn_2015 <- glPca(gl_2015_syn_rmNA, useC = TRUE, nf=20)

save(pca_syn_1990, file="pca_syn_1990.Rdata")
save(pca_syn_2015, file="pca_syn_2015.Rdata")

#Load for subsequent analysis if need be
#load("pca_syn_1990.Rdata")
#load("pca_syn_2015.Rdata")

#Perform DAPC in adegenet
#Set pop by BT, MN-BT (where we know for sure there is BT), and NON-BT, rather than by just states, since these are our biologically informative "groups"
setPop(gl_1990_syn_rmNA) <- ~BT
setPop(gl_2015_syn_rmNA) <- ~BT

#Next, run dapc with all PCs retained to examine a-score at each level. 
dapc_syn_1990 <- dapc(gl_1990_syn_rmNA, gl_1990_syn_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_syn_1990)
temp_syn_1990 <- optim.a.score(dapc_syn_1990) #Use 9 PCs 
dapc_syn_1990 <- dapc(gl_1990_syn_rmNA, gl_1990_syn_rmNA@pop, n.pca = 9, n.da = 3, glPca = pca_syn_1990)

dapc_syn_2015 <- dapc(gl_2015_syn_rmNA, gl_2015_syn_rmNA@pop, n.pca = 20, n.da = 3, glPca = pca_syn_2015)
temp_syn_2015 <- optim.a.score(dapc_syn_2015) #Use 11 PCs
dapc_syn_2015 <- dapc(gl_2015_syn_rmNA, gl_2015_syn_rmNA@pop, n.pca = 11, n.da = 3, glPca = pca_syn_2015)

#Make plots of results
##1990

#Make a function to save plots
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

myCol = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")
par(bty = "n") #Remove box around plots
scatter(dapc_syn_1990, lwd = 2, lty = 2, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = myCol, xlim = c(-4,4))
points(dapc_syn_1990$grp.coord[,1], dapc_syn_1990$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_syn_1990$grp.coord[,1], dapc_syn_1990$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=myCol)
scatter_syn_1990 <- grab_grob()

#Reformat so that individuals are sorted by their expected membership so that compoplot is easier to interpret
compoplot_syn_1990_df <- data.frame(grp=dapc_syn_1990$grp, dapc_syn_1990$posterior)
compoplot_syn_1990_df <- compoplot_syn_1990_df[order(compoplot_syn_1990_df$grp),]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_syn_1990_df[,2:4]), cex.axis = 0.75, cex.names = 0.5, col.pal = myCol, show.lab = TRUE, legend = FALSE)
compo_syn_1990 <- grab_grob()

##2015
scatter(dapc_syn_2015, lwd = 2, lty = 2, scree.da = FALSE, legend = FALSE, clabel = FALSE, cstar=0, solid=.6, col = c("BT" = "#D55E00", "NO BT" = "#E69F00", "MN BT" = "#56B4E9"))
legend(3, 2.5, c("NO BT", "MN BT", "BT"), bty = "n", pch = 20, col = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00"))
points(dapc_syn_2015$grp.coord[,1], dapc_syn_2015$grp.coord[,2], pch=4, cex=1.5, lwd=3, col="black")
points(dapc_syn_2015$grp.coord[,1], dapc_syn_2015$grp.coord[,2], pch=4, cex=1.5, lwd=1, col=c("BT" = "#D55E00", "NO BT" = "#E69F00", "MN BT" = "#56B4E9"))
scatter_syn_2015 <- grab_grob()

compoplot_syn_2015_df <- data.frame(grp=dapc_syn_2015$grp, dapc_syn_2015$posterior)
#Re-order factor levels and columns to be same order as 1990
compoplot_syn_2015_df$grp = factor(compoplot_syn_2015_df$grp,levels(compoplot_syn_2015_df$grp)[c(2,3,1)])
compoplot_syn_2015_df <- compoplot_syn_2015_df[order(compoplot_syn_2015_df$grp),]
compoplot_syn_2015_df <- compoplot_syn_2015_df[,c(1,3,4,2)]
par(mai=c(1.5,1,1,1))
compoplot(as.matrix(compoplot_syn_2015_df[,2:4]), col.pal = myCol, cex.axis = 0.75, cex.names = 0.5, show.lab = TRUE, legend = FALSE)
compo_syn_2015 <- grab_grob()

#Save the plots
pdf("1990_2015_syn_dapc.pdf", height = 9, width = 7)
grid.arrange(arrangeGrob(scatter_syn_1990, top = "1990", clip = "on"), arrangeGrob(scatter_syn_2015, top = "2015", clip = "on"), arrangeGrob(compo_syn_1990, clip = "on"), arrangeGrob(compo_syn_2015, clip = "on"), ncol=2, nrow=2, vp = viewport(width=0.98, height=0.98))
dev.off()




####Running analysis with k=2 or 3 to compare results#########
#Rather than using already defined biological groups, one can also identify clusters de novo
#Plots can be examined for "elbow" in the chart to identify the correct k, and below is an example of this with our data
#grp_1990 <- find.clusters(gl_1990_rmNA, stat = "BIC", glPca = pca_1990, n.pca = 10, max.n.clust=20) #elbow at k=13
#grp_2015 <- find.clusters(gl_2015_rmNA, stat = "BIC", glPca = pca_2015,  n.pca = 9, max.n.clust=20) #elbow at k=12
#grp_1990_2015 <- find.clusters(gl_1990_2015_rmNA, stat = "BIC", glPca = pca_1990_2015, n.pca = 15, max.n.clust=20) #elbow at k=17
#Then, DAPC can be performed with these groups rather than biological groups
#dapc_1990k <- dapc(gl_1990_rmNA, grp_1990$grp, n.pca = 10, n.da = 3, glPca = pca_1990)
#dapc_2015k <- dapc(gl_2015_rmNA, grp_2015$grp, n.pca = 9, n.da = 3, glPca = pca_2015)
#dapc_1990_2015k <- dapc(gl_1990_2015_rmNA, grp_1990_2015$grp, n.pca = 15, n.da = 3, glPca = pca_1990_2015)

#Here, we are testing with just k=2 or 3, depending on dataset, to see how to it compares to that number of biological groups
grp_1990 <- find.clusters(gl_1990_rmNA, stat = "BIC", glPca = pca_1990, n.pca = 20, max.n.clust=20) #chose k = 3
dapc_1990k3 <- dapc(gl_1990_rmNA, grp_1990$grp, n.pca = 20, n.da = 3, glPca = pca_1990)
temp_1990 <- optim.a.score(dapc_1990k3) #Use 3 PCs since the spline interpolation suggested that as the ideal amount.
dapc_1990k3 <- dapc(gl_1990_rmNA, grp_1990$grp, n.pca = 3, n.da = 3, glPca = pca_1990)

grp_2015 <- find.clusters(gl_2015_rmNA, stat = "BIC", glPca = pca_2015, n.pca = 20, max.n.clust=20) #chose k = 3
dapc_2015k3 <- dapc(gl_2015_rmNA, grp_2015$grp, n.pca = 20, n.da = 3, glPca = pca_2015)
temp_2015 <- optim.a.score(dapc_2015k3) #Use 4 PCs since the spline interpolation suggested that as the ideal amount.
dapc_2015k3 <- dapc(gl_2015_rmNA, grp_2015$grp, n.pca = 4, n.da = 3, glPca = pca_2015)

grp_1990_2015 <- find.clusters(gl_1990_2015_rmNA, stat = "BIC", glPca = pca_1990_2015, n.pca = 20, max.n.clust=20) #chose k = 2
dapc_1990_2015k2 <- dapc(gl_1990_2015_rmNA, grp_1990_2015$grp, n.pca = 20, n.da = 3, glPca = pca_1990_2015)
temp_1990_2015 <- optim.a.score(dapc_1990_2015k2) #Use 1 PCs since the spline interpolation suggested that as the ideal amount.
dapc_1990_2015k2 <- dapc(gl_1990_2015_rmNA, grp_1990_2015$grp, n.pca = 1, n.da = 3, glPca = pca_1990_2015)

#Make tables of membership probabilities for k = 2 and k = 3 vs defined based on geography

#Reformat so that individuals are sorted by their expected membership so that compoplot is easier to interpret
#These plots show membership probabilities of each individual for the different biological groups based on the retained discriminant functions
#Membership probabilities also provide indications of how clear-cut genetic clusters are.
#Loose clusters will result in fairly flat distributions of membership probabilities of individuals across clusters, pointing to possible admixture.
compoplot_1990k3_df <- data.frame("Group" = dapc_1990k3$grp, "Membership Probabilities Group 1" = dapc_1990k3$posterior[,1], "Membership Probabilities Group 2" = dapc_1990k3$posterior[,2], "Membership Probabilities Group 3" = dapc_1990k3$posterior[,3])
compoplot_1990k3_df <- compoplot_1990k3_df[order(compoplot_1990k3_df$Group),]

compoplot_2015k3_df <- data.frame("Group" = dapc_2015k3$grp, "Membership Probabilities Group 1" = dapc_2015k3$posterior[,1], "Membership Probabilities Group 2" = dapc_2015k3$posterior[,2], "Membership Probabilities Group 3" = dapc_2015k3$posterior[,3])
compoplot_2015k3_df <- compoplot_2015k3_df[order(compoplot_2015k3_df$Group),]

compoplot_1990_2015k2_df <- data.frame("Group" = dapc_1990_2015k2$grp, "Membership Probabilities Group 1" = dapc_1990_2015k2$posterior[,1], "Membership Probabilities Group 2" = dapc_1990_2015k2$posterior[,2])
compoplot_1990_2015k2_df <- compoplot_1990_2015k2_df[order(compoplot_1990_2015k2_df$Group),]

library("openxlsx")
write.xlsx(rbind(compoplot_1990k3_df, compoplot_2015k3_df), "k3_years_separate.xlsx", rowNames = TRUE)
#Then, add population information in to Excel file

write.xlsx(compoplot_1990_2015k2_df, "k2_years_together.xlsx", rowNames = TRUE)
