## Make the test vector for controlling for latitude and longitude in hgdp data set
## This script takes two sets of .afreq files and subsets each to variants over 1% and then gets the intersection of the two lists

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript make_Tvec_hgdp_cordinates.R <fam file> <populations> <samples> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

fam_file = args[1]
pop_file = args[2]
sample_file = args[3]
out_file = args[4]
