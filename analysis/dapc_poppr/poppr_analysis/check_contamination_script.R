#This script is for data analysis with vcfR using the information in the recent publication (https://www.frontiersin.org/articles/10.3389/fgene.2018.00123/full) to examine isolate contamination

#Change working directory to appropriate location
setwd("")

library("vcfR") #1.8.0
library("ggplot2")

#Read in VCF files for 1990 and 2015 isolates and save as Rdata
vcf_1990 <- read.vcfR("../adegenet_dapc_analysis/1990_isolates.filter.vcf")
vcf_2015 <- read.vcfR("../adegenet_dapc_analysis/2015_isolates.filter.vcf")
#vcf_1990_2015 <- read.vcfR("../adegenet_dapc_analysis/1990_2015_isolates.filter.vcf")

save(vcf_1990, file="vcf_1990.Rdata")
save(vcf_2015, file="vcf_2015.Rdata")
#save(vcf_1990_2015, file="vcf_1990_2015.Rdata")

#Then reload objects for subsequent analyses if need be
load("vcf_1990.Rdata")
load("vcf_2015.Rdata")
#load("vcf_1990_2015.Rdata")

#Subset to het positions
#This subsetting will place an NA in non-heterozygous positions in the vcfR object
vcf_1990_gt <- extract.gt(vcf_1990, element = 'GT')
vcf_1990_hets <- is_het(vcf_1990_gt)
is.na(vcf_1990@gt[,-1][ !vcf_1990_hets ]) <- TRUE

vcf_2015_gt <- extract.gt(vcf_2015, element = 'GT')
vcf_2015_hets <- is_het(vcf_2015_gt)
is.na(vcf_2015@gt[,-1][ !vcf_2015_hets ]) <- TRUE

#Filter het positions on allele depth
vcf_1990_ad <- extract.gt(vcf_1990, element = 'AD')
vcf_2015_ad <- extract.gt(vcf_2015, element = 'AD')

vcf_1990_allele1 <- masplit(vcf_1990_ad, record = 1)
vcf_1990_allele2 <- masplit(vcf_1990_ad, record = 2)

vcf_2015_allele1 <- masplit(vcf_2015_ad, record = 1)
vcf_2015_allele2 <- masplit(vcf_2015_ad, record = 2)

#Produce quantiles at 0.15 and 0.95 probabilites
sums_1990 <- apply(vcf_1990_allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)
sums_2015 <- apply(vcf_2015_allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)

#Allele 1
vcf_1990_dp2 <- sweep(vcf_1990_allele1, MARGIN=2, FUN = "-", sums_1990[1,])
vcf_2015_dp2 <- sweep(vcf_2015_allele1, MARGIN=2, FUN = "-", sums_2015[1,])

#allele1[dp2 < 0] <- NA
vcf_1990@gt[,-1][ vcf_1990_dp2 < 0 & !is.na(vcf_1990@gt[,-1]) ] <- NA
vcf_1990_dp2 <- sweep(vcf_1990_allele1, MARGIN=2, FUN = "-", sums_1990[2,])

vcf_2015@gt[,-1][ vcf_2015_dp2 < 0 & !is.na(vcf_2015@gt[,-1]) ] <- NA
vcf_2015_dp2 <- sweep(vcf_2015_allele1, MARGIN=2, FUN = "-", sums_2015[2,])

#allele1[dp2 > 0] <- NA
vcf_1990@gt[,-1][vcf_1990_dp2 > 0] <- NA
vcf_2015@gt[,-1][vcf_2015_dp2 > 0] <- NA

#Allele 2
vcf_1990_dp2 <- sweep(vcf_1990_allele2, MARGIN=2, FUN = "-", sums_1990[1,])
vcf_1990@gt[,-1][ vcf_1990_dp2 < 0 & !is.na(vcf_1990@gt[,-1]) ] <- NA
vcf_1990_dp2 <- sweep(vcf_1990_allele2, MARGIN=2, FUN = "-", sums_1990[2,])
vcf_1990@gt[,-1][vcf_1990_dp2 > 0] <- NA

vcf_2015_dp2 <- sweep(vcf_2015_allele2, MARGIN=2, FUN = "-", sums_2015[1,])
vcf_2015@gt[,-1][ vcf_2015_dp2 < 0 & !is.na(vcf_2015@gt[,-1]) ] <- NA
vcf_2015_dp2 <- sweep(vcf_2015_allele2, MARGIN=2, FUN = "-", sums_2015[2,])
vcf_2015@gt[,-1][vcf_2015_dp2 > 0] <- NA

#Calculate allele balance with the new filtered vcf
vcf_1990_ad <- extract.gt(vcf_1990, element = 'AD')
vcf_2015_ad <- extract.gt(vcf_2015, element = 'AD')

vcf_1990_allele1 <- masplit(vcf_1990_ad, record = 1)
vcf_1990_allele2 <- masplit(vcf_1990_ad, record = 2)

vcf_2015_allele1 <- masplit(vcf_2015_ad, record = 1)
vcf_2015_allele2 <- masplit(vcf_2015_ad, record = 2)

vcf_1990_ad1 <- vcf_1990_allele1 / (vcf_1990_allele1 + vcf_1990_allele2)
vcf_1990_ad2 <- vcf_1990_allele2 / (vcf_1990_allele1 + vcf_1990_allele2)

vcf_2015_ad1 <- vcf_2015_allele1 / (vcf_2015_allele1 + vcf_2015_allele2)
vcf_2015_ad2 <- vcf_2015_allele2 / (vcf_2015_allele1 + vcf_2015_allele2)

#Make a tidy dataframe to make plot
vcf_1990_ad1t <- tidyr::gather(tibble::as.tibble(vcf_1990_ad1), "Sample", "Ab")
vcf_1990_ad1t$Allele <- "ab1"

vcf_1990_ad2t <- tidyr::gather(tibble::as.tibble(vcf_1990_ad2), "Sample", "Ab")
vcf_1990_ad2t$Allele <- "ab2"

vcf_1990_t <- rbind(vcf_1990_ad1t, vcf_1990_ad2t)
vcf_1990_t <- dplyr::filter(vcf_1990_t, !is.na(Ab))

#Remove reference isolates for 2015 before plotting
vcf_2015_ad1t <- tidyr::gather(tibble::as.tibble(vcf_2015_ad1), "Sample", "Ab")
vcf_2015_ad1t <- vcf_2015_ad1t[!(vcf_2015_ad1t$Sample %in% c("12SD80", "12NC29")), ]
vcf_2015_ad1t$Allele <- "ab1"

vcf_2015_ad2t <- tidyr::gather(tibble::as.tibble(vcf_2015_ad2), "Sample", "Ab")
vcf_2015_ad2t <- vcf_2015_ad2t[!(vcf_2015_ad2t$Sample %in% c("12SD80", "12NC29")), ]
vcf_2015_ad2t$Allele <- "ab2"

vcf_2015_t <- rbind(vcf_2015_ad1t, vcf_2015_ad2t)
vcf_2015_t <- dplyr::filter(vcf_2015_t, !is.na(Ab))

#Make plot
p_1990 <- ggplot(vcf_1990_t, aes(x=Ab))
p_1990 <- p_1990 + scale_x_continuous(breaks=c(0, 1/2, 1), labels = c('0', '1/2', '1'))
p_1990 <- p_1990 + geom_histogram(data=subset(vcf_1990_t, Allele == 'ab1'),fill =  "#1f78b4", binwidth = 0.02)
p_1990 <- p_1990 + geom_histogram(data=subset(vcf_1990_t, Allele == 'ab2'),fill =  "#a6cee3", binwidth = 0.02)
p_1990 <- p_1990 + facet_wrap( ~ Sample, ncol=10, nrow = 3)
p_1990 <- p_1990 + theme_bw()

p_1990 <- p_1990 + theme(panel.grid.major.x=element_line(color = "#A9A9A9", size=0.4) )
p_1990 <- p_1990 + theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.1, linetype = 'dashed') )
p_1990 <- p_1990 + theme(panel.grid.minor.x=element_blank())
p_1990 <- p_1990 + theme(panel.grid.minor.y=element_blank())
p_1990 <- p_1990 + theme(strip.text.x = element_text(size = 8))
p_1990 <- p_1990 + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8))
p_1990 <- p_1990 + theme(axis.text.y = element_text(angle = 30, hjust = 1, size=10))

p_1990 <- p_1990 + theme(axis.title.x = element_text(size = 14, angle = 0))
p_1990 <- p_1990 + theme(axis.title.y = element_text(size = 14, angle = 90))
p_1990 <- p_1990 + xlab("Allele balance")
p_1990 <- p_1990 + ylab("Heterozygous positions")

#Save the plot
pdf("allele_balance_1990.pdf", height = 9, width = 9)
p_1990
dev.off()

#Make plot
p_2015 <- ggplot(vcf_2015_t, aes(x=Ab))
p_2015 <- p_2015 + scale_x_continuous(breaks=c(0, 1/2, 1), labels = c('0', '1/2', '1'))
p_2015 <- p_2015 + geom_histogram(data=subset(vcf_2015_t, Allele == 'ab1'),fill =  "#1f78b4", binwidth = 0.02)
p_2015 <- p_2015 + geom_histogram(data=subset(vcf_2015_t, Allele == 'ab2'),fill =  "#a6cee3", binwidth = 0.02)
p_2015 <- p_2015 + facet_wrap( ~ Sample, ncol=10, nrow = 3)
p_2015 <- p_2015 + theme_bw()

p_2015 <- p_2015 + theme(panel.grid.major.x=element_line(color = "#A9A9A9", size=0.4) )
p_2015 <- p_2015 + theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.1, linetype = 'dashed') )
p_2015 <- p_2015 + theme(panel.grid.minor.x=element_blank())
p_2015 <- p_2015 + theme(panel.grid.minor.y=element_blank())
p_2015 <- p_2015 + theme(strip.text.x = element_text(size = 8))
p_2015 <- p_2015 + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8))
p_2015 <- p_2015 + theme(axis.text.y = element_text(angle = 30, hjust = 1, size=10))

p_2015 <- p_2015 + theme(axis.title.x = element_text(size = 14, angle = 0))
p_2015 <- p_2015 + theme(axis.title.y = element_text(size = 14, angle = 90))
p_2015 <- p_2015 + xlab("Allele balance")
p_2015 <- p_2015 + ylab("Heterozygous positions")

#Save the plot
pdf("allele_balance_2015.pdf", height = 9, width = 9)
p_2015
dev.off()
