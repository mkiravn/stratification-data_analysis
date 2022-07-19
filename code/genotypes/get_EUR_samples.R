## Get EUR samples
## This script subsets the list of HGDP individuals to only those in europe

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript get_EUR_sample.R <population meta data> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

popfile = args[1]
outfile = args[2]


# Read in metadata
pops <- fread(popfile)

# Subset to Europeans
eur <- pops %>% filter(region == "EUROPE")

# Format output
eur <- eur %>% select(sample)
eur$SEX <- rep(NA, nrow(eur))
colnames(eur) <- c("#IID", "SEX")

# Write output
fwrite(eur, outfile,row.names = F, col.names = T, quote = F, sep = "\t")
