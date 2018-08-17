#This script is for analysis of r squared LD values calculated by vcftools
#I got inspiration from these two posts for this script:
#https://www.biostars.org/p/300381/
#https://www.biostars.org/p/84443/

library(dplyr)
library(stringr)
library(ggplot2)

#Change working directory to appropriate location
setwd("")

r2_1990_BT <- read.delim("1990_isolates.LD.200kb.BT.geno.ld", header = TRUE, colClasses = c("character", rep("numeric", 4)))
r2_1990_NOBT <- read.delim("1990_isolates.LD.200kb.NOBT.geno.ld", header = TRUE, colClasses = c("character", rep("numeric", 4)))
r2_1990_MNBT <- read.delim("1990_isolates.LD.200kb.MNBT.geno.ld", header = TRUE, colClasses = c("character", rep("numeric", 4)))

r2_2015_BT <- read.delim("2015_isolates.LD.200kb.BT.geno.ld", header = TRUE, colClasses = c("character", rep("numeric", 4)))
r2_2015_NOBT <- read.delim("2015_isolates.LD.200kb.NOBT.geno.ld", header = TRUE, colClasses = c("character", rep("numeric", 4)))
r2_2015_MNBT <- read.delim("2015_isolates.LD.200kb.MNBT.geno.ld", header = TRUE, colClasses = c("character", rep("numeric", 4)))

#Calculate the distance between snp1 and snp2 in kb
r2_1990_BT$distancekb <- (r2_1990_BT$POS2-r2_1990_BT$POS1)/1000
r2_1990_NOBT$distancekb <- (r2_1990_NOBT$POS2-r2_1990_NOBT$POS1)/1000
r2_1990_MNBT$distancekb <- (r2_1990_MNBT$POS2-r2_1990_MNBT$POS1)/1000

r2_2015_BT$distancekb <- (r2_2015_BT$POS2-r2_2015_BT$POS1)/1000
r2_2015_NOBT$distancekb <- (r2_2015_NOBT$POS2-r2_2015_NOBT$POS1)/1000
r2_2015_MNBT$distancekb <- (r2_2015_MNBT$POS2-r2_2015_MNBT$POS1)/1000

#Bin 200 kb
r2_1990_BT$grp <- cut(r2_1990_BT$distancekb, 0:200)
r2_1990_NOBT$grp <- cut(r2_1990_NOBT$distancekb, 0:200)
r2_1990_MNBT$grp <- cut(r2_1990_MNBT$distancekb, 0:200)

r2_2015_BT$grp <- cut(r2_2015_BT$distancekb, 0:200)
r2_2015_NOBT$grp <- cut(r2_2015_NOBT$distancekb, 0:200)
r2_2015_MNBT$grp <- cut(r2_2015_MNBT$distancekb, 0:200)

##r2 mean every 1kb
r2_1990_means_BT <- r2_1990_BT %>% group_by(grp) %>% summarise(mean = mean(R.2, na.rm = TRUE), median = median(R.2, na.rm = TRUE))
r2_1990_means_NOBT <- r2_1990_NOBT %>% group_by(grp) %>% summarise(mean = mean(R.2, na.rm = TRUE), median = median(R.2, na.rm = TRUE))
r2_1990_means_MNBT <- r2_1990_MNBT %>% group_by(grp) %>% summarise(mean = mean(R.2, na.rm = TRUE), median = median(R.2, na.rm = TRUE))

r2_2015_means_BT <- r2_2015_BT %>% group_by(grp) %>% summarise(mean = mean(R.2, na.rm = TRUE), median = median(R.2, na.rm = TRUE))
r2_2015_means_NOBT <- r2_2015_NOBT %>% group_by(grp) %>% summarise(mean = mean(R.2, na.rm = TRUE), median = median(R.2, na.rm = TRUE))
r2_2015_means_MNBT <- r2_2015_MNBT %>% group_by(grp) %>% summarise(mean = mean(R.2, na.rm = TRUE), median = median(R.2, na.rm = TRUE))

#A helper step to get mid points of distance intervals for plotting
r2_1990_means_BT_plot <- r2_1990_means_BT %>% mutate(start = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")), end = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")), mid=start+((end-start)/2))

r2_1990_means_NOBT_plot <- r2_1990_means_NOBT %>% mutate(start = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")), end = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")), mid=start+((end-start)/2))

r2_1990_means_MNBT_plot <- r2_1990_means_MNBT %>% mutate(start = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")), end = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")), mid=start+((end-start)/2))

r2_2015_means_BT_plot <- r2_2015_means_BT %>% mutate(start = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")), end = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")), mid=start+((end-start)/2))

r2_2015_means_NOBT_plot <- r2_2015_means_NOBT %>% mutate(start = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")), end = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")), mid=start+((end-start)/2))

r2_2015_means_MNBT_plot <- r2_2015_means_MNBT %>% mutate(start = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")), end = as.integer(str_extract(str_replace_all(grp,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")), mid=start+((end-start)/2))

#Set theme to be classic
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=12),
              axis.title = element_text(size=15),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 13)
            )
)

#Make plots
ld_1990 <- ggplot() +
  geom_point(data=r2_1990_means_BT_plot, aes(x=start,y=mean), size=0.4, colour="#D55E00") +
  geom_line(data=r2_1990_means_BT_plot, aes(x=start,y=mean), size=0.3, alpha=0.5,colour="#D55E00") +
  geom_point(data=r2_1990_means_NOBT_plot, aes(x=start,y=mean), size=0.4, colour="#E69F00") +
  geom_line(data=r2_1990_means_NOBT_plot, aes(x=start,y=mean), size=0.3, alpha=0.5,colour="#E69F00") +
  geom_point(data=r2_1990_means_MNBT_plot, aes(x=start,y=mean), size=0.4, colour="#56B4E9") +
  geom_line(data=r2_1990_means_MNBT_plot, aes(x=start,y=mean), size=0.3, alpha=0.5,colour="#56B4E9") +
  labs(x="Distance (Kilobases)", y=expression(LD~(r^{2}))) +
  scale_x_continuous(breaks=c(0,50,100,150,200), labels=c("0","50","100","150","200")) +
  geom_hline(yintercept=0.2, color = "red", size=0.5)

ld_2015 <- ggplot() +
  geom_point(data=r2_2015_means_BT_plot, aes(x=start,y=mean), size=0.4, colour="#D55E00") +
  geom_line(data=r2_2015_means_BT_plot, aes(x=start,y=mean), size=0.3, alpha=0.5,colour="#D55E00") +
  geom_point(data=r2_2015_means_NOBT_plot, aes(x=start,y=mean), size=0.4, colour="#E69F00") +
  geom_line(data=r2_2015_means_NOBT_plot, aes(x=start,y=mean), size=0.3, alpha=0.5,colour="#E69F00") +
  geom_point(data=r2_2015_means_MNBT_plot, aes(x=start,y=mean), size=0.4, colour="#56B4E9") +
  geom_line(data=r2_2015_means_MNBT_plot, aes(x=start,y=mean), size=0.3, alpha=0.5,colour="#56B4E9") +
  labs(x="Distance (Kilobases)", y=expression(LD~(r^{2}))) +
  scale_x_continuous(breaks=c(0,50,100,150,200), labels=c("0","50","100","150","200")) +
  geom_hline(yintercept=0.2, color = "red", size=0.5)

pdf("r2_1990.pdf")
ld_1990
dev.off()

pdf("r2_2015.pdf")
ld_2015
dev.off()
