#In this script, the GC content of shared and unique variants between different callers is explored

library(stringi)
library(dplyr)
library(ggplot2)
library(ggridges)

#Set working directory to appropriate location
setwd("~/Documents/Work/Log Book Data/2016/2016_crownrust_NIFA_project/sequencing_data_analysis/gatk_freebayes_samtools_compare/gc_analysis")

#Place all files with .cov extension in vector
#Make sure to use double backslash to escape the period in the filename
#Only read in 2015 and 1990 files, make sure to exclude any 2012 in the directory
gc_files <- list.files(pattern="*.gc.*$", full.names=TRUE)

##Read in files and rename columns
#Rename columns vector
gc_cols <- c("contig", "start", "end", "NULL", "Percent GC", rep("NULL", 7))
#Make a function to read in files
loadGCFile <- function(x) {
  #read in relevant columns from a gc file, extract the sample name from the file, and add it as a column
  #also add in a column if 1990, 2012, or 2015 isolate
  df <- read.delim(x, header=TRUE, colClasses = c("character", "integer", "integer", "NULL", "numeric", rep("NULL", 7)), col.names=gc_cols)
  # df$id <- sub("[21590]{4}_gatk_freebayes_samtools_isolates.+", "Shared", basename(x))
  # df$id <- sub("[21590]{4}_isolates_unique_to_freebayesVSall.+", "Unique FreeBayes", basename(x))
  # df$id <- sub("[21590]{4}_isolates_unique_to_gatkVSall.+", "Unique GATK", basename(x))
  # df$id <- sub("[21590]{4}_isolates_unique_to_samtoolsVSall.+", "Unique Samtools", basename(x))
  df$id <- stri_replace_all_regex(str = basename(x), pattern = c("[21590]{4}_gatk_freebayes_samtools_isolates.+", "[21590]{4}_isolates_unique_to_freebayesVSall.+", "[21590]{4}_isolates_unique_to_gatkVSall.+", "[21590]{4}_isolates_unique_to_samtoolsVSall.+"), replacement = c("Shared", "Unique FreeBayes", "Unique GATK", "Unique Samtools"), vectorize_all = FALSE)
  df$Year <- as.factor(sub("([21590]{4}).+", "\\1", basename(x)))
  df$Variant <- ifelse(grepl("snp", basename(x)), "SNP", "INDEL")
  df
}

#Then use lapply to apply this function to all files and load in as list
gc <- lapply(gc_files, loadGCFile)

#Then, make this list into one long dataframe
gc_df <- bind_rows(gc)

#Dot plots
#Set theme to be classic
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=12),
              axis.title = element_text(size=15),
              axis.line = element_line(size=0),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 13)
            )
)

#Ridge plots
plot <- ggplot(gc_df, aes(x=(gc_df$Percent.GC)*100, y=id)) +
               geom_density_ridges(aes(fill = Year), color = "white", alpha = 0.5, scale=1) + 
               scale_fill_manual(values=c("black", "red")) +
               scale_y_discrete(expand=c(0.001, 0)) +
               scale_x_continuous(limits = c(20, 75)) +
               facet_wrap(~Variant) +
               labs(x = "Percent GC Surrounding Variants", y="")

#Save the plot
pdf("variant_gc_neighborhoods.pdf")
plot
dev.off()
