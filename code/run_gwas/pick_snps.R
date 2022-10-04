## Pick SNPs
## This script reads in .glm.linear files (uncorrected, lat, long), assigns SNPs to LD blocks (pickrell 2016) and picks the SNP with the lowest p-value per block

args=commandArgs(TRUE)

if(length(args)<8){stop("Rscript pick_snps.R <blocks> <uncorrected betas> <lat corrected betas> <long corrected betas> <outfile uncorrected> <outfile lat> <outfile long>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

block_file = args[1]
u_file = args[2]
lat_file = args[3]
long_file = args[4]
out_u = args[5]
out_lat = args[6]
out_long = args[7]
pt = as.numeric(args[8])
print(pt)

# Function to assign snps to ld blocks
assign_SNP_to_block <- function(betas) {

  # Add new column for LD block
  betas$block <- NA

  # Get chromosome of interest
  block_chr <- ld_blocks %>% filter(chr == as.numeric(betas$`#CHROM`[1]))

  # Get first SNP actually in LD block (skips SNPs out of range at the begining)
  if (betas$POS[1] < as.numeric(block_chr[1,2])) {
    i <- tail(which(betas$POS < as.numeric(block_chr[1,2])),1)+1
  } else {
    i <- 1
  }
  bp <- betas$POS[i]

  # Loop through LD blocks
  for (block_num in block_chr$block_number) {
    block_num <- as.numeric(as.character(block_num))
    s <- as.numeric(block_chr %>% filter(block_number == block_num) %>% select(start))
    e <- as.numeric(block_chr %>% filter(block_number == block_num) %>% select(stop))
    betas <- betas %>% mutate(block = case_when((POS >= s & POS < e) ~ block_num, TRUE ~ as.numeric(as.character(block))))
  }
  return(betas)
}

# Read in block file
ld_blocks = fread(block_file)

# Read in all betas and convert to correct type
betas_u <- fread(u_file)
betas_u$P <- as.numeric(betas_u$P)
betas_lat <- fread(lat_file)
betas_lat$P <- as.numeric(betas_lat$P)
betas_long <- fread(long_file)
betas_long$P <- as.numeric(betas_long$P)
print("Read in data")

# Eliminate SNPs above p-value threshold
df_u <- betas_u %>%
  filter(P < pt)
df_lat <- betas_lat %>%
  filter(P < pt)
df_long <- betas_long %>%
  filter(P < pt)
print("filtered SNPs")

# Assign remaining SNPs to block
df_u <- assign_SNP_to_block(df_u)
df_lat <- assign_SNP_to_block(df_lat)
df_long <- assign_SNP_to_block(df_long)

# Pick the minimum p-value per block
u <- df_u %>%
  drop_na() %>% group_by(block) %>% arrange(P) %>% slice(n=1)
lat <- df_lat %>%
  drop_na() %>% group_by(block) %>% arrange(P) %>% slice(n=1)
long <- df_long %>%
  drop_na() %>% group_by(block) %>% arrange(P) %>% slice(n=1)

# Write to output files
fwrite(u,out_u, row.names = F, col.names = T, quote = F, sep = "\t")
print("wrote uncorrected")
fwrite(lat,out_lat, row.names = F, col.names = T, quote = F, sep = "\t")
print("wrote latidue")
fwrite(long,out_long, row.names = F, col.names = T, quote = F, sep = "\t")
print("wrote longitude")


