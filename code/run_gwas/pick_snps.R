## Pick SNPs
## This script reads in .glm.linear files (uncorrected, lat, long), assigns SNPs to LD blocks (pickrell 2016) and picks the SNP with the lowest p-value per block

args=commandArgs(TRUE)

if(length(args)<8){stop("Rscript pick_snps.R <blocks> <uncorrected betas> <lat corrected betas> <long corrected betas> <outfile uncorrected> <outfile lat> <outfile long>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  #library(tidyr)
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

# Read in block file
ld_blocks = fread(block_file)


# Function to assign SNPs to blocl
assign_SNP_to_block <- function(CHR, BP, block = ld_blocks) {

  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start)

  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}

# Read in all betas
betas_u <- fread(u_file)
betas_lat <- fread(lat_file)
betas_long <- fread(long_file)


test <- betas_u %>% filter(P < pt) 

# Threshold and asign all SNPs to blocks
df_u <- betas_u %>%
  filter(P < pt) %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[1]), as.numeric(params[2]))))
df_lat <- betas_lat %>%
  filter(P < pt) %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[1]), as.numeric(params[2]))))
df_long <- betas_long %>%
  filter(P < pt) %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[1]), as.numeric(params[2]))))


# Pick the minimum p-value per block
u <- df_u %>%
  group_by(block) %>% arrange(P) %>% slice(n=1)
lat <- df_lat %>%
  group_by(block) %>% arrange(P) %>% slice(n=1)
long <- df_long %>%
  group_by(block) %>% arrange(P) %>% slice(n=1)


# Write to output files
fwrite(u,out_u, row.names = F, col.names = T, quote = F, sep = "\t")
print("wrote uncorrected")
fwrite(lat,out_lat, row.names = F, col.names = T, quote = F, sep = "\t")
print("wrote latidue")
fwrite(long,out_long, row.names = F, col.names = T, quote = F, sep = "\t")
print("wrote longitude")


