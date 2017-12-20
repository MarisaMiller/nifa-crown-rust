#This R script was substantially modified from the visualiztion R script that is downloaded with kWIP (see https://kwip.readthedocs.io/en/latest/)
#!/usr/bin/env Rscript
library(fields, quietly = T)
setwd("")

dist.file = paste0("2015", ".dist")
kern.file = paste0("2015", ".kern")

dist.file = paste0("1990", ".dist")
kern.file = paste0("1990", ".kern")

dist.file = paste0("1990_2015", ".dist")
kern.file = paste0("1990_2015", ".kern")

base = basename("2015")
base = basename("1990")
base = basename("1990_2015")

dist = as.matrix(read.delim(dist.file)[, -1])

#2015
#Don't need to include 2012 isolates
dist <- dist[-c(2:3), -c(2:3)]

colnames(dist) <- c("15FL1.2", "15Fl1.4", "15MN10.4", "15MN10.5", "15MN13.4", "15MN14.4", "15MN15.3", "15MN16.3", "15MN17.5", "15MN18.1", "15MN18.3", "15MN23.1", "15MN24.1", "15MN25.3", "15MN27.3", "15MS7.1", "15ND19.2", "15ND19.5", "15ND20.3", "15ND20.4", "15NE8.4", "15NE8.5", "15NE9.1", "15NE9.3", "15OH12.3", "15SD11.1", "15SD11.2", "15SD30.1", "15SD30.3", "15TX3.1")

#1990
#colnames(dist) <- c("90AR100-1",	"90GA16-1",	"90KS101-1",	"90LA38-1",	"90MN137-1",	"90MN13B-3",	"90MN148-1",	"90MN149-2",	"90MN14B-1",	"90MN152-1",	"90MN153-1",	"90MN17B-1",	"90MN1B-1",	"90MN2B-1",	"90MN3B-1",	"90MN5B-1",	"90MN7B-1",	"90MN8B-2",	"90MN9B-4",	"90PA162-1",	"90SD164-1",	"90SD171-1",	"90SD172-1",	"90TX45-1",	"90TX47-1",	"90TX52-1",	"90TX58-1",	"90TX70-1",	"90WI131-1",	"90WI132-1")

#1990_2015
#Don't need to include 2012 isolates
#dist <- dist[-c(2:3), -c(2:3)]

#colnames(dist) <- c("15FL1.2", "15Fl1.4", "15MN10.4", "15MN10.5", "15MN13.4", "15MN14.4", "15MN15.3", "15MN16.3", "15MN17.5", "15MN18.1", "15MN18.3", "15MN23.1", "15MN24.1", "15MN25.3", "15MN27.3", "15MS7.1", "15ND19.2", "15ND19.5", "15ND20.3", "15ND20.4", "15NE8.4", "15NE8.5", "15NE9.1", "15NE9.3", "15OH12.3", "15SD11.1", "15SD11.2", "15SD30.1", "15SD30.3", "15TX3.1", "90AR100-1",	"90GA16-1",	"90KS101-1",	"90LA38-1",	"90MN137-1",	"90MN13B-3",	"90MN148-1",	"90MN149-2",	"90MN14B-1",	"90MN152-1",	"90MN153-1",	"90MN17B-1",	"90MN1B-1",	"90MN2B-1",	"90MN3B-1",	"90MN5B-1",	"90MN7B-1",	"90MN8B-2",	"90MN9B-4",	"90PA162-1",	"90SD164-1",	"90SD171-1",	"90SD172-1",	"90TX45-1",	"90TX47-1",	"90TX52-1",	"90TX58-1",	"90TX70-1",	"90WI131-1",	"90WI132-1")

kern = as.matrix(read.delim(kern.file)[, -1])

#Don't need to include 2012 isolates for 2015 and combined analysis
kern <- kern[-c(2:3), -c(2:3)]

n = ncol(dist)

clust = hclust(as.dist(dist))
sam_ord = clust$order

#Below are ways to plot from the original R script with kWIP
# pdf(paste0(base, ".pdf"))
# image.plot(1:n, 1:n, dist[sam_ord, sam_ord],
#            axes=TRUE,
#            xlab="",
#            ylab="",
#            col=rainbow(80, start=0.2),
#            main=paste("Distance of", base))
# axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.5, las=2, tick=F)
# axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.5, las=2, tick=F)
# 
# image.plot(1:n, 1:n, kern[sam_ord, sam_ord],
#            axes=FALSE,
#            xlab="",
#            ylab="",
#            col=rainbow(80, start=0.2),
#            main=paste("Kernel of", base))
# axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.5, las=2, tick=F)
# axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.5, las=2, tick=F)
# 
# d = diag(kern)
# geom.mean.d = sqrt(d %*% t(d))
# normk = kern / geom.mean.d
# 
# image.plot(1:n, 1:n, normk[sam_ord, sam_ord],
#            axes=FALSE,
#            xlab="",
#            ylab="",
#            col=rainbow(80, start=0.2),
#            main=paste("Normalised Kernel of", base))
# axis(1, 1:n, colnames(dist)[sam_ord], srt=45, xpd=T, cex.axis=0.5, las=2, tick=F)
# axis(2, 1:n, colnames(dist)[sam_ord], srt=90, cex.axis=0.5, las=2, tick=F)
# 
# plot(clust, cex=0.5, main=base)

fit <- cmdscale(dist, eig=TRUE, k=2) # k is the number of dim

# plot solution
# x <- fit$points[,1]
# y <- fit$points[,2]
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
#   main="Metric MDS", type='p')
# text(x, y, labels = colnames(dist), cex=.4)
# 
# dev.off()

#Use ggplot
#First, merge by rows
library(plyr)

a <- as.data.frame(fit$points[,1])
b <- as.data.frame(fit$points[,2])
names <- as.data.frame(colnames(dist))

a$rn <- rownames(a)
b$rn <- rownames(b)
names$rn <- rownames(names)

df <- join_all(list(a, b, names), by = 'rn', type = 'full')

library(ggplot2)

#Add rownames back
rownames(df) <- c(as.character(df$`colnames(dist)`))

#2015
df$IsolateState <- as.factor(c("FL", "FL", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MS", "ND", "ND", "ND", "ND", "NE", "NE", "NE",  "NE", "OH", "SD", "SD", "SD", "SD", "TX"))

#1990
#df$IsolateState <- as.factor(c("AR",	"GA",	"KS",	"LA",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"PA",	"SD",	"SD",	"SD",	"TX",	"TX",	"TX",	"TX",	"TX",	"WI",	"WI"))

#1990_2015
#df$IsolateState <- as.factor(c("FL", "FL", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MN", "MS", "ND", "ND", "ND", "ND", "NE", "NE", "NE",  "NE", "OH", "SD", "SD", "SD", "SD", "TX", "AR",	"GA",	"KS",	"LA",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"MN",	"PA",	"SD",	"SD",	"SD",	"TX",	"TX",	"TX",	"TX",	"TX",	"WI",	"WI"))

#df$IsolateYear <- as.factor(c("2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015", "2015",  "2015", "2015", "2015", "2015", "2015", "2015", "2015", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990", "1990",  "1990", "1990", "1990", "1990", "1990", "1990", "1990"))

pdf("kWIP_2015.pdf")
#pdf("kWIP_1990.pdf")
#pdf("kWIP_1990_2015.pdf")

ggplot(df, aes(x=df$`fit$points[, 1]`, y=df$`fit$points[, 2]`, color = df$IsolateState)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("MN" = "firebrick2", "SD" = "gold", "TX" = "dodgerblue", "ND" = "black", "NE" = "darkorchid", "OH" = "darkgreen", "MS" = "darkorange", "FL" = "brown", "WI" = "darkgrey", "PA" = "deeppink", "KS" = "darkolivegreen4", "AR" = "darkorange3", "GA" = "tomato3", "LA" = "slateblue2")) +
  theme_classic(base_size = 16) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab("Coordinate 1") +
  ylab("Coordinate 2")

#For 1990_2015
# ggplot(df, aes(x=df$`fit$points[, 1]`, y=df$`fit$points[, 2]`, color = df$IsolateYear)) +
#   geom_point(size = 2) +
#   scale_color_manual(values = c("2015" = "firebrick2", "1990" = "black")) +
#   theme_classic(base_size = 16) +
#   theme(legend.title=element_blank(), legend.position="top") +
#   xlab("Coordinate 1") +
#   ylab("Coordinate 2")

dev.off()

#Older plot options
#library(ggrepel)
#ggplot(df) +
#  geom_point(aes(df$`fit$points[, 1]`, df$`fit$points[, 2]`), color = "grey", size = 3) +
#geom_label_repel(aes(df$`fit$points[, 1]`, df$`fit$points[, 2]`, fill = IsolateState, label = rownames(df)),
#                 fontface = 'bold', color = 'white',
#                 box.padding = 0.35, point.padding = 0.5,
#                 segment.color = 'grey50') +
#  theme_classic(base_size = 16) +
#  scale_fill_manual(guide = guide_legend(override.aes = aes(color = NA)), values = setNames(c("brown", "firebrick2", "darkorange", "azure4", "black", "darkorchid", "darkgreen", "gold", "dodgerblue"), levels(df$IsolateState))) +
#  xlab("Coordinate 1") +
#  ylab("Coordinate 2") +
#  labs(fill="Isolate State")