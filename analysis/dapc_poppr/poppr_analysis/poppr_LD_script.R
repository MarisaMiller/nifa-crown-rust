#This script is for data analysis with poppr using the information at https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html

#Change working directory to appropriate location
setwd("")

#Only need vcfR for reading in original vcf files and doing some filtering
library("vcfR") #1.8.0
library("poppr") #2.7.1
library("ggplot2")
#library("ggridges") optional, if you decide to use ridgeline plots instead of boxplots

#These are the same steps that were taken for the DAPC analysis in dapc_script.R in the adegenet_dapc_analysis directory
#Read in VCF files for 1990 and 2015 isolates, just do once and then convert to genlight (see below)
#vcf_1990 <- read.vcfR("1990_isolates.filter.vcf")
#vcf_2015 <- read.vcfR("2015_isolates.filter.vcf")

#Convert the large and unwieldy VCF files into smaller genlight objects that can be used with poppr/adegenet
#Non-biallelic variants will be discarded since genlight objects don't support non-biallelic data
#gl_1990 <- vcfR2genlight(vcf_1990)
#gl_2015 <- vcfR2genlight(vcf_2015)

#save(gl_1990, file="gl_1990.Rdata")
#save(gl_2015, file="gl_2015.Rdata")

#Then reload objects for subsequent analyses if need be
# load("gl_1990.Rdata")
# load("gl_2015.Rdata")

#Remove reference isolates (12NC29 and 12SD80) for this analysis and drop alleles no longer polymorphic in the data
#gl_2015 <- gl_2015[!(indNames(gl_2015) %in% c("12SD80", "12NC29"))]

#Set population information to samples, note that this is different from the DAPC analysis script
#pop(gl_1990) <- as.factor(c("NO BT", "NO BT", "NO BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "MN BT", "MN BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "MN BT", "BT", "BT", "MN BT", "NO BT", "BT", "BT", "NO BT", "NO BT", "BT", "NO BT", "MN BT"))

#pop(gl_2015) <- as.factor(c("BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "BT", "BT", "BT", "BT", "NO BT", "BT", "NO BT", "BT", "MN BT", "BT", "MN BT", "BT", "MN BT", "BT", "NO BT", "BT", "MN BT"))

#Set ploidy
#gl_list <- list(gl_1990, gl_2015)
#lapply(gl_list, ploidy, 2)

#First remove any NAs
# toRemove1990 <- is.na(glMean(gl_1990, alleleAsUnit = FALSE)) # TRUE where NA
# which(toRemove1990) # position of entirely non-typed loci
# gl_1990_rmNA_LD <- gl_1990[, !toRemove1990]
# toRemove2015 <- is.na(glMean(gl_2015, alleleAsUnit = FALSE))
# which(toRemove2015)
# gl_2015_rmNA_LD <- gl_2015[, !toRemove2015]

# save(gl_1990_rmNA_LD, file="gl_1990_rmNA_LD.Rdata")
# save(gl_2015_rmNA_LD, file="gl_2015_rmNA_LD.Rdata")

#Then reload objects for subsequent analyses if need be
load("../adegenet_dapc_analysis/gl_1990_rmNA_LD.Rdata") #1,379,292 biallelic SNPs used for analysis
load("../adegenet_dapc_analysis/gl_2015_rmNA_LD.Rdata") #1,186,114 biallelic SNPs

#First look at the index of association in randomly sampled sites (10,000 SNP sets repeated 100 times)
#This function will return the standardized index of association ("rbarD")
ia_1990 <- samp.ia(gl_1990_rmNA_LD, n.snp = 10000L, reps = 100L)
ia_2015 <- samp.ia(gl_2015_rmNA_LD, n.snp = 10000L, reps = 100L)

#Means
#1990
#0.061
#2015
#0.053

#Generate simulated data with different levels of clonality to test against for each year
#Performed as in https://apsjournals.apsnet.org/doi/full/10.1094/MPMI-10-17-0258-R
#See code I used for inspiration at https://github.com/grunwaldlab/rubi-gbs/blob/master/Step5_Data_processing.Rmd

#In simulations using 30 individuals and 1,282,703 biallelic SNP loci (average between 1990 and 2015 loci amounts)
### No structure (fully admixed pops)
sex <- glSim(30, 1282703, ploid=2, LD=T)
### Semi-clonal 
clone_50 <- glSim(30, 1282703, n.snp.struc=641352, ploid=2, LD=T)
### Most-clonal 
clone_75 <- glSim(30, 1282703, n.snp.struc=962027, ploid=2, LD=T)
### Structure (clonal pops)
clone_100 <- glSim(30, 1282703, n.snp.struc=1282703, ploid=2, LD = T)

## IA sex
ia.sex <- samp.ia(sex, n.snp = 10000L, reps = 100L)
## IA clone
ia.clone.50 <- samp.ia(clone_50, n.snp = 10000L, reps = 100L)
## IA.semiclone
ia.clone.75 <- samp.ia(clone_75, n.snp = 10000L, reps = 100L)
## IA.mostclone
ia.clone.100 <- samp.ia(clone_100, n.snp = 10000L, reps = 100L)

# Summarizing data frames
d1 <- data.frame(ia_1990, rep("1990_dataset", length(ia_1990)))
d2 <- data.frame(ia_2015, rep("2015_dataset", length(ia_2015)))
d3 <- data.frame(ia.sex, rep("sexual", length(ia.sex)))
d4 <- data.frame(ia.clone.50, rep("clone_50", length(ia.clone.50)))
d5 <- data.frame(ia.clone.75, rep("clone_75", length(ia.clone.75)))
d6 <- data.frame(ia.clone.100, rep("clone_100", length(ia.clone.100)))
colnames(d1) <- c("ia","dset")
colnames(d2) <- c("ia","dset")
colnames(d3) <- c("ia","dset")
colnames(d4) <- c("ia","dset")
colnames(d5) <- c("ia","dset")
colnames(d6) <- c("ia","dset")
ia.total <- rbind(d1, d2, d3, d4, d5, d6)

# Normality tests
#The Shapiro-Wilk’s normality test showed that the following observed and simulated distributions did not follow normality (P < 0.05): 1990_ia, 2015_ia, 50% simulated linkage dataset, 75% simulated linkage dataset, and the 100% linkage set. The other set (0%) did show normality.
frames <- list(as.data.frame(d1), as.data.frame(d2), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5), as.data.frame(d6))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

#A nonparametric multiple comparison of ranks using Kruskal-Wallis was performed to test for differences across mean ranks.
#The P. rubi sample mean was situated between the simulated 75% linkage data and the simulated 50% linkage data. This indicated that populations of P. rubi are predominantly clonal and/or selfing.
library(agricolae)

# Kruskal wallis test
kruskal.test(ia ~ dset, ia.total) #Kruskal-Wallis chi-squared = 577.07, df = 5, p-value < 2.2e-16
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))
# ia groups
# clone_100    546.41      a
# 1990_dataset 454.59      b
# 2015_dataset 346.29      c
# clone_75     254.71      d
# clone_50     150.50      e
# sexual        50.50      f

#Split by pop in each year and repeat analysis
no_bt_1990 <- popsub(gl_1990_rmNA_LD, "NO BT")
bt_1990 <- popsub(gl_1990_rmNA_LD, "BT")
mn_bt_1990 <- popsub(gl_1990_rmNA_LD, "MN BT")

no_bt_2015 <- popsub(gl_2015_rmNA_LD, "NO BT")
bt_2015 <- popsub(gl_2015_rmNA_LD, "BT")
mn_bt_2015 <- popsub(gl_2015_rmNA_LD, "MN BT")

ia_no_bt_1990 <- samp.ia(no_bt_1990, n.snp = 10000L, reps = 100L)
ia_bt_1990 <- samp.ia(bt_1990, n.snp = 10000L, reps = 100L)
ia_mn_bt_1990 <- samp.ia(mn_bt_1990, n.snp = 10000L, reps = 100L)

ia_no_bt_2015 <- samp.ia(no_bt_2015, n.snp = 10000L, reps = 100L)
ia_bt_2015 <- samp.ia(bt_2015, n.snp = 10000L, reps = 100L)
ia_mn_bt_2015 <- samp.ia(mn_bt_2015, n.snp = 10000L, reps = 100L)

#Means
#1990
#NOBT BT  MNBT
#0.146  0.035  0.112
#2015
#NOBT BT  MNBT
#0.269  0.063 0.005

#Make dataframe for plotting and stats
ia_no_bt_1990_df <- data.frame(as.data.frame(ia_no_bt_1990), "NO BT", 1990)
colnames(ia_no_bt_1990_df) <- c("rbarD", "Population", "Year")

ia_bt_1990_df <- data.frame(as.data.frame(ia_bt_1990), "BT", 1990)
colnames(ia_bt_1990_df) <- c("rbarD", "Population", "Year")

ia_mn_bt_1990_df <- data.frame(as.data.frame(ia_mn_bt_1990), "MN BT", 1990)
colnames(ia_mn_bt_1990_df) <- c("rbarD", "Population", "Year")

ia_no_bt_2015_df <- data.frame(as.data.frame(ia_no_bt_2015), "NO BT", 2015)
colnames(ia_no_bt_2015_df) <- c("rbarD", "Population", "Year")

ia_bt_2015_df <- data.frame(as.data.frame(ia_bt_2015), "BT", 2015)
colnames(ia_bt_2015_df) <- c("rbarD", "Population", "Year")

ia_mn_bt_2015_df <- data.frame(as.data.frame(ia_mn_bt_2015), "MN BT", 2015)
colnames(ia_mn_bt_2015_df) <- c("rbarD", "Population", "Year")

ia_pop_df <- rbind(ia_no_bt_1990_df, ia_bt_1990_df, ia_mn_bt_1990_df, ia_no_bt_2015_df, ia_bt_2015_df, ia_mn_bt_2015_df)
ia_pop_df_1990 <- rbind(ia_no_bt_1990_df, ia_bt_1990_df, ia_mn_bt_1990_df)
ia_pop_df_2015 <- rbind(ia_no_bt_2015_df, ia_bt_2015_df, ia_mn_bt_2015_df)

#Kruskal wallis test
kruskal.test(rbarD ~ Population, ia_pop_df) #Kruskal-Wallis chi-squared = 399.33, df = 2, p-value < 2.2e-16

k.test <- with(ia_pop_df_1990, kruskal(rbarD, Population, group = T, p.adj = "bon"))
# $groups
# rbarD groups
# NO BT 250.5      a
# MN BT 150.5      b
# BT     50.5      c

k.test <- with(ia_pop_df_2015, kruskal(rbarD, Population, group = T, p.adj = "bon"))
# $groups
# rbarD groups
# NO BT 250.5      a
# BT    150.5      b
# MN BT  50.5      c


#Plot distributions for each year and populations
#Set theme to be classic for rest of plots
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=10),
              axis.title = element_text(size=15),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 13)
            )
)

byYear_plot <- ggplot(ia.total, aes(x = dset, y = ia, color = as.factor(dset))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("1990", "2015", "Sexual", "50% Linkage", "75% Linkage", "100% Linkage")) +
  scale_y_continuous(limits = c(0, 0.075)) +
  scale_color_manual(values = c("black", "red", rep("grey", 4))) +
  theme(legend.position="none") +
  labs(x = "", y = expression(italic(bar("r")["d"])))

# byYear_plot <- ggplot(ia_year_df, aes(x = rbarD, y = Year, fill = as.factor(Year))) +
#   #stat_density_ridges(quantile_lines = TRUE) +
#   geom_density_ridges() +
#   scale_fill_manual(values = c("black", "red")) +
#   theme(legend.title = element_blank()) +
#   scale_y_discrete(expand = c(0.001, 0)) +
#   labs(x = expression(italic(bar("r")["d"])), y = "")

#Save the plots
pdf("ia_byYear.pdf")
byYear_plot
dev.off()

byPop_plot <- ggplot(ia_pop_df, aes(x = Population, y = rbarD, color = as.factor(Population))) +
  geom_boxplot() +
  scale_color_manual(values = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")) +
  theme(legend.title = element_blank(), panel.grid.major.y = element_line(color = "black"), axis.text.x = element_blank()) +
  labs(x = "Population", y = expression(italic(bar("r")["d"]))) +
  facet_wrap(~Year)

# byPop_plot <- ggplot(ia_pop_df, aes(x = rbarD, y = Population, fill = as.factor(Population))) +
#   #stat_density_ridges(quantile_lines = TRUE) +
#   scale_color_manual(values = c("black")) +
#   geom_density_ridges(scale = 4) +
#   scale_fill_manual(values = c("NO BT" = "#E69F00", "MN BT" = "#56B4E9", "BT" = "#D55E00")) +
#   theme(legend.title = element_blank(), panel.grid.major.y = element_line(color = "black")) +
#   scale_y_discrete(expand = c(0.001, 0)) +
#   labs(x = expression(italic(bar("r")["d"])), y = "Population") +
#   facet_wrap(~Year)

#Save the plots
pdf("ia_byPop.pdf")
byPop_plot
dev.off()


####Sliding window based analysis####
#Check with new version first
# remotes::install_github("grunwaldlab/poppr@positiion")
# library(poppr)
# 
# win_ia_1990 <- win.ia(gl_1990_rmNA_LD, window = 10000L, min.snps = 3L)
# win_ia_2015 <- win.ia(gl_2015_rmNA_LD, window = 10000L, min.snps = 3L)
