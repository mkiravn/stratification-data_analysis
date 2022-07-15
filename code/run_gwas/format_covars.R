## Format covar
## This script reads in all the potential covariates (age at reccruitment, genetic sex, array type, TGWAS) and formats them for plink

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript conccat_Tm.R <prefix to Tm chromosomes> <chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

TGWAS = fread(args[1])
aar = fread(args[2])
sex = fread(args[3])
array = fread(args[4])
outfile = args[5]
