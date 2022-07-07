## Overlapping SNPs
## This script takes two sets of .afreq files and subsets each to variants over 1% and then gets the intersection of the two lists

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript overlapping_snps.R <ukbb .freq> <test panel .freq>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

ukbb_file = args[1]
tp_file = args[2]

# Load freuqncy files
ukbb <- fread(ukbb_file)
tp <- fread(tp_file)

# Subset to SNPs > 1% MAF
ukbb <- subset(ukbb, ukbb$ALT_FREQS > 0.01 & ukbb$ALT_FREQS < 0.99)
tp <- subset(tp, tp$ALT_FREQS > 0.01 & tp$ALT_FREQS < 0.99)

#
test <- inner_join(ukbb, tp, by = c("ID", "ALT", "REF"))

