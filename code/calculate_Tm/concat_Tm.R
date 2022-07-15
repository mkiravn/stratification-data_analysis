## Concatenate projections
## This script reads in all the projected test vectors for each chromosome and adds them all together and scales the final covariate

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript conccat_Tm.R <prefix to Tm chromosomes> <chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
outfile = args[2]

chrs <- rep(NA, (length(args) - 2))
for (j in 1:length(chrs)) {
    chrs[j] <- as.numeric(args[2+j])
}
chrs <- sort(chrs)

# Add Tm for each chromosome to each other
df <- fread(paste0(Tm_prefix, "_", chrs[1], ".txt"))
for (i in 2:length(chrs)) {

  new <- fread(paste0(Tm_prefix, "_", chrs[i], ".txt"))
  df$latitude <- df$latitude + new$latitude
  df$longitude <- df$longitude + new$longitude

}

df$latitude <- scale(df$latitude)
df$longitude <- scale(df$longitude)

# Save output
fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







