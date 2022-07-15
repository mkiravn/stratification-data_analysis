## Concatenate projections
## This script reads in all the projected test vectors for each chromosome and adds them all together and scales the final covariate

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript conccat_Tm.R <prefix to Tm chromosomes> <chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
chr = args[2]
outfile = args[3]

print(head(chr))


# Add Tm for each chromosome to each other
df <- fread(paste0(Tm_prefix, "_22.txt"))
for (i in 21:1) {

  new <- fread(paste0(Tm_prefix, "_", i, ".txt"))
  df$latitude <- df$latitude + new$latitude
  df$longitude <- df$longitude + new$longitude

}

# Save output
fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







