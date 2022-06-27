# Script formats standing height phenotype file

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide path to phenotype file")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

pheno_file = args[1]
out_file = args[2]

# Read in phenotypes
df <- fread(pheno_file)
df <-fread("~/stratification-data_analysis/standing_height_50.txt")

# Select ID and initial visit heights and rename columns
df <- df %>% select(`50-0.0`, `50-1.0`)
colnames(df) <- c("IID", "Height")

# Remove individuals with no initial value for height
df <- df %>% drop_na(Height)

# Save file
fwrite(df, out_file,col.names=T,row.names=F,quote=F,sep="\t")
