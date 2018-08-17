#This script is for calculating genome coverage in bins for a batch of files post read-mapping
#Change to the appropriate working directory where all coverage files are stored

#Load appropriate libraries
library(dplyr)
library(ggplot2)

#Place all files with .cov extension in vector
#Make sure to use double backslash to escape the period in the filename
cov_files <- list.files(pattern="^[19].*\\.cov$", full.names=TRUE)

##Read in files and rename columns
#Rename columns vector
covcols <- c("contig", "start", "end", "contig_window", "depth", "number of bases at depth", "window_size", "percent of window at depth")

#Make a function to read in files
loadCovFile <- function(x) {
  #read in a coverage file, extract the sample name from the file, and add it as a column
  #also add in a column if 1990 or 2015 isolate
  df <- read.delim(x, header=FALSE, col.names=covcols)
  df$sample_name <- sub("(.+)\\.window\\.cov", "\\1", basename(x))
  df$sample_year <- sub("([12590]{2}).+\\.cov", "\\1", basename(x))
  df
}

#Then use lapply to apply this function to all files and load in as list
cov <- lapply(cov_files, loadCovFile)

#Then, make this list into one long dataframe
cov_df <- do.call(rbind, cov)

#Remove "all" rows as they're not useful for this analysis
cov_df <- subset(cov_df, contig != "all")

#Filter out windows that have more than 75% of the window at a depth of less than 2.
#Then, sum up the length of 50kb windows that are "missing" to estimate PAV.
#Of course this is influenced by sequencing depth, so it's hard to determine if these parts of genome are truly missing or not.
missing_windows <- cov_df %>% group_by(.dots=c("sample_name", "contig_window")) %>%
  filter(depth < 2) %>%
  summarize(fraction_window_missing = sum(percent.of.window.at.depth)) %>%
  filter(fraction_window_missing > 0.75) %>%
  count(sample_name) %>%
  mutate(total_missing_windows_kb = (n*50)) %>%
  arrange(desc(total_missing_windows_kb))
  
#Format for saving as table so that coverage of individual windows can be queried
missing_windows_comprehensive <- cov_df %>% group_by(.dots=c("sample_name", "contig_window", "start", "end")) %>%
  filter(depth < 2) %>%
  summarize(fraction_window_missing = sum(percent.of.window.at.depth))

write.table(missing_windows_comprehensive, file = "fraction_missing_per_window.txt", sep = "\t", row.names = FALSE)

#Can also do on a per-base basis
missing_bases <- cov_df %>% filter(depth < 2) %>%
  group_by(sample_name) %>%
  summarize(missing_bases_Mb = sum(number.of.bases.at.depth)/1000000) %>%
  arrange(desc(missing_bases_Mb))


#MN5B is interesting! It has higher than average homozygosity, and has a fairly unique virulence pattern.
#Does it have a large deletion in the genome?

#now make plot
#Set theme to be classic
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=10),
              axis.title = element_text(size=15),
              axis.line = element_line(size=0.5),
              strip.text.y = element_blank()
            )
)

#Make plot
plot_cov <- ggplot(missing_windows, aes(x=reorder(sample_name, -total_missing_windows_kb), y=total_missing_windows_kb)) +
                geom_point() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(x = "", y="Sum of 50 kb windows with low depth")

# Save the plots
pdf("window_coverage_plot.pdf", width = 8, height = 8)
plot_cov
dev.off()
