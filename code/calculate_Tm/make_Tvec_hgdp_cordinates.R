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
out_file_ALL = args[3]
out_file_EUR = args[4]

# Read in files
fam <- fread(fam_file)
pops <- fread(pop_file)

# Merge sample and fam files
df <- inner_join(fam, pops, by = c("#IID"= "sample"))
df <- df %>% select("#IID", "latitude", "longitude")

# Mean center coordinates
df$latitude <- df$latitude - mean(df$latitude)
df$longitude <- df$longitude - mean(df$longitude)

# Output file
fwrite(df,paste0(out_file_ALL, "_cordinates.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# Output separate phenotype file for lat and long
lat <- df %>% select("#IID", "latitude")
fwrite(lat,paste0(out_file_ALL, "_latitude.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

long <- df %>% select("#IID", "longitude")
fwrite(long,paste0(out_file_ALL, "_longitude.txt"), row.names = F, col.names = T, quote = F, sep = "\t")


##### EUR only

# Merge sample and fam files
df <- inner_join(fam, pops, by = c("#IID"= "sample"))
df <- df %>% filter(region == "EUROPE") %>% select("#IID", "latitude", "longitude")

# Mean center coordinates
df$latitude <- df$latitude - mean(df$latitude)
df$longitude <- df$longitude - mean(df$longitude)

# Output file
fwrite(df,paste0(out_file_EUR, "_cordinates.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# Output separate phenotype file for lat and long
lat <- df %>% select("#IID", "latitude")
fwrite(lat,paste0(out_file_EUR, "_latitude.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

long <- df %>% select("#IID", "longitude")
fwrite(long,paste0(out_file_EUR, "_longitude.txt"), row.names = F, col.names = T, quote = F, sep = "\t")





