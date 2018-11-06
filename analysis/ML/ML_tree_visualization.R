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
