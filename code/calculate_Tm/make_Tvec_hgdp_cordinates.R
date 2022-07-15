## Make the test vector for controlling for latitude and longitude in hgdp data set
## This script takes two sets of .afreq files and subsets each to variants over 1% and then gets the intersection of the two lists

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript make_Tvec_hgdp_cordinates.R <fam file> <populations> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

fam_file = args[1]
pop_file = args[2]
out_file = args[3]

# Read in files
fam <- fread(fam_file)
pops <- fread(pop_file)

# Merge sample and fam files
df <- inner_join(fam, pops, by = c("#IID"= "sample"))
df <- df %>% select("#IID", "SEX", "latitude", "longitude")

# Mean center coordinates
df$latitude <- df$latitude - mean(df$latitude)
df$longitude <- df$longitude - mean(df$longitude)

# Output file
fwrite(df,out_file, row.names = F, col.names = T, quote = F, sep = "\t")










