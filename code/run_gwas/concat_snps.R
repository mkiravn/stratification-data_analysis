## concat_snps.R
## This script reads in ascertained SNPs for all chromosomes separately and concatenates them into the same file

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript concat_snps.R <prefix for ascertained snps> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

snp_prefix = args[1]
outfile_u = args[2]
outfile_lat = args[3]
outfile_long = args[4]

read_betas <- function(beta_suffix) {

  # Read in betas
  df <- fread(paste0(snp_prefix, 1, beta_suffix))

  for (i in 2:22) {

    # Read in betas
    betas <- fread(paste0(snp_prefix, i, beta_suffix))

    # Add to previous df
    df <- rbind(df,  betas)

  }

  # Flip betas to get the effect size of the ALT allele
  df <- df %>% mutate(BETA_Strat = case_when(ALT == A1 ~ BETA, REF == A1 ~ -1 * BETA))

  # Return Betas
  return(df)
}

# Concat Betas
df_u <- read_betas("_v3.Height.betas")
df_lat <- read_betas("_v3.Height-Lat.betas")
df_long <- read_betas("_v3.Height-Long.betas")

# Write files
fwrite(df_u, outfile_u,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_lat, outfile_lat,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_long, outfile_long,row.names=F,quote=F,sep="\t", col.names = T)


