#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
################################################## CROWN RUST PHYLOGENY ###############################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# Maximum Likelihood Phylogenetics for Crown Rust 1990 and 2015 data
# Melaria Figueroa, Marisa Miller, Hoa Nguyen-Phuc

###################################################### DATA SET ##################################################

# 61 samples: 60 samples in 1990 and 2015 PLUS only one 2012 sample (12NC29)
# VCF file is 1990_2015_isolates.filter.12NC29.vcf

# This file was prepared as follows:
# You can create a PBS script or just run quickly at the command line
VCFLIB=/home/kianians/millerme/nifa-crown-rust/scripts/vcflib/bin
module load freebayes/20161103 parallel liblzma/5.2.2 gcc/7.2.0

cd $PBS_O_WORKDIR

$VCFLIB/vcfremovesamples ../snp_calling/freebayes/1990_2015_isolates.filter.vcf 12SD80 | $VCFLIB/vcffixup - > 1990_2015_isolates.filter.12NC29.vcf

#######################################################################################################################
########################################################### PART 1 ####################################################
################################################## VCF TO DNAbin TO PHYLIP ############################################
#######################################################################################################################

# THIS IS R CODE

# Use R packages "vcfR" to convert VCF to DNAbin and then "ips" (RAxML wrapper) to convert DNAbin to PHYLIP

# Convert VCF to DNAbin ----------------------------------------------------------------------------------------------

library(vcfR) #1.8.0

# Read VCF file
mainvcf = read.vcfR("1990_2015_isolates.filter.12NC29.vcf") 

# Convert VCF to DNAbin
maindna = vcfR2DNAbin(mainvcf,
                      extract.haps = FALSE,
                      unphased_as_NA = FALSE,
                      #ref.seq = ref, # this is optional, turn it on with small dataset
                      consensus = TRUE,
                      verbose = TRUE)

# Convert DNAbin to PHYLIP -------------------------------------------------------------------------------------------

library(ips)
raxmlpath = ("D:\\GoogleDrive\\400PseudoCodes\\420Packages\\raxml\\raxmlHPC") ## need to specify executable, the path on your machine may be different!

MLtree = try(raxml(maindna, 
                   f = "a", m="GTRCAT", N = "autoMRE", 
                   p = 1234, x = 1234, threads = 4, 
                   file= "DNAbin", exec = raxmlpath))

# This will: 
# 1) generate 2 ALIGNED FASTA files - "DNAbin.phy"and "DNAbin.phy.reduced" -  specified in file "RAxML_info.DNAbin":
# 2) run 10,000 bootstraps (default) 
# the .reduced is the one with undetermined heterozygote sites removed. Use this file!!

# The run will take a really long time (about 40m for 1 bootstrap) using personal PC machine
# stop the run (only takes a few minutes to make .PHY files)
# take the .PHY file to run in your supercomputing Linux platform


#######################################################################################################################
########################################################### PART 2 ####################################################
############################################## RAxML MAXIMUM LIKELIHOOD INFERENCES ####################################
#######################################################################################################################

# Create PBS file - save as "raxml.pbs"
#!/bin/bash -l
#PBS -l walltime=69:59:00,nodes=1:ppn=24,mem=60gb
#PBS -m abe 
#PBS -M hnguyenp@umn.edu
cd $PBS_O_WORKDIR
module load raxml/8.2.11_pthread
raxml -s DNAbin.phy.reduced -n data2.no.outgroup -f a -m GTRCAT -p 12345 -x 12345 -# 500 --no-bfgs -T 24

#######################################################################################################################
######################################################### PART 3 ######################################################
################################################## RAXML TOUCHING-UP ##################################################
#######################################################################################################################

# To test if there are enough bootstraps

module load raxml/8.2.11_pthread
raxml -m GTRCAT -z RAxML_bootstrap.data2.no.outgroup -I autoMRE -p 12345 -n data2.no.outgroup 

#######################################################################################################################
######################################################### PART 4 ######################################################
################################################ MAKE PHYLOGENTIC TREES ###############################################
#######################################################################################################################

# This is R code

# Read "best" ML tree with support Bootstrap values and branch (edge) length records:
# Note: change the RAxML file here accordingly, it should start with "RAxML_bipartitions"
# and end with your pbs submitting name, 
# e.g. job with "rax3.1ref.pbs" will have output name "RAxML_bipartitions.rax3.1ref"

library("ggtree")
library("ape")
library("phangorn")
library("ggplot2")

BestTreeML = read.tree("RAxML_bestTree.data2.no.outgroup")
AnnoTreeML = read.tree("RAxML_bipartitions.data2.no.outgroup")
RootTreeML = midpoint(AnnoTreeML)    
newick = RootTreeML

#Make a dataframe with color information for isolates from populations based on alternate host status
no_bt <- c("90GA16-1","90TX52-1","90TX47-1","90KS101-1","90TX45-1","90LA38-1","90TX58-1","90AR100-1","90TX70-1", "15MS7-1","15TX3-1","15FL1-2","15Fl1-4")
bt <- c("90MN137-1","90MN148-1","90MN149-2","90MN152-1","90MN153-1","90WI131-1","90WI132-1","90SD164-1","90SD171-1","90SD172-1","90PA162-1", "15ND19-2","15ND19-5","15ND20-3","15ND20-4","15MN10-4","15MN10-5","15MN27-3","15MN23-1","15MN24-1","15MN25-3","15SD30-1","15SD30-3","15SD11-1","15SD11-2","15NE8-4","15NE8-5","15NE9-1","15NE9-3","15OH12-3")
mn_bt <- c("90MN1B-1","90MN2B-1","90MN3B-1","90MN5B-1","90MN7B-1","90MN8B-2","90MN9B-4","90MN13B-3","90MN14B-1","90MN17B-1", "15MN13-4","15MN14-4","15MN15-3","15MN16-3","15MN17-5","15MN18-1","15MN18-3")

popColorInfo <- data.frame(isolate = newick$tip.label, population = NA, color = NA)
popColorInfo$population[popColorInfo$isolate %in% no_bt] <- "NO BT"
popColorInfo$population[popColorInfo$isolate %in% bt] <- "BT"
popColorInfo$population[popColorInfo$isolate %in% mn_bt] <- "MN BT"

popColorInfo$color[popColorInfo$population == "NO BT"] <- "#E69F00"
popColorInfo$color[popColorInfo$population == "BT"] <- "#D55E00"
popColorInfo$color[popColorInfo$population == "MN BT"] <- "#56B4E9"

groupInfo <- split(newick$tip.label, gsub("([12590]{2}).+", "\\1", newick$tip.label))
newick <- groupOTU(newick, groupInfo)

tree_plot <- ggtree(newick, aes(color = group), layout='circular') + geom_tiplab(size=3, aes(angle=angle), color = "black", offset = 0.01) + geom_point2(aes(label=label, subset=!isTip & as.numeric(label) < 50), color = "green", size=1) 

tree_plot <- tree_plot %<+% popColorInfo + geom_tippoint(aes(color = color), size=1) + scale_color_manual(values = c("#56B4E9", "#D55E00", "#E69F00", "purple", "red", "black")) + geom_treescale(x=.35, y=0)

tree_plot <- flip(tree_plot, 87, 83)

#Save plot
pdf("ml_tree.pdf", height = 7, width = 7)
tree_plot
dev.off()

#######################################################################################################################
#######################################################################################################################
########################################################### END #######################################################
#######################################################################################################################
#######################################################################################################################
