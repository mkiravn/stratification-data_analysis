## Pick SNPs
## This script reads in ascertained SNPs and the test vector and computes Qx and an empirical p-value

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript compute Qx.R <prefix for ascertained snps> <test panel prefix> <test vec file> <qx outfile> <pgs outfile>   ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

snp_prefix = args[1]
genos_prefix = args[2]
tvec_file = args[3]
outfile_qx = args[4]
outfile_pgs = args[5]

chr_start=21
chr_end=22


# Function to read in genotype matrix for a set of variants (counted allele is the Alternate)
read_genos <- function(geno_prefix, betas) {

  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas$ID
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=F)
  colnames(X) <- var.ids

  return(X)
}

# Function to calculate Va
calc_Va <- function(geno_mat, es) {

  # Get allele frequency in test panel
  freq <- colMeans(geno_mat) /2

  # Pull out effect sizes only
  effect_size <- es$BETA_Strat

  # Compute Va
  Va <- 2 * sum((effect_size)^2 * freq * (1 - freq))
  return(Va)
}


# Function to compute PGS
pgs <- function(X, betas) {

  # Mean center genotypes
  X <- apply(X, 2, function(y) y - mean(y))

  # Comput PGS
  Bhat_strat <- betas$BETA_Strat
  Z_strat <- X %*% Bhat_strat

  # Format output
  out <- cbind(Z_strat)
  colnames(out) <- c("STRAT")

  return(out)
}

# Function to calculate Qx
calc_Qx <- function(mprs, tvec, Va, lambda_T) {

  # Compute Qx Strat
  Ztest <- t(tvec) %*% mprs$strat.adjusted
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va*lambda_T)

  return(Qx_strat)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(betas, tvec, Va, X, true_file, lambda_T) {

  # Flip effect sizes
  betas$BETA_Strat <- flip(betas$BETA_Strat)

  # Calculate PGS
  prs <- pgs(X, betas)

  # Calculate Qx
  Qx <- t(calc_Qx(stand_PGS(prs, true_file), tvec, Va, lambda_T))

  return(Qx)
}



# Load in Test vectors
TV <- fread(tvec_file)








# Calcluate PGS by looping through chromosomes
N <- nrow(fread(paste0(genos_prefix, i, ".psam")))
vecVa <- rep(0, (chr_end - chr_start + 1))
matPGS <- matrix(0,nrow = N,  ncol = (chr_end - chr_start + 1))
for (i in chr_start:chr_end) {

  # Read in betas
  betas <- fread(paste0(snp_prefix, i))

  # Flip betas to get the effect size of the ALT allele
  betas <- betas %>% mutate(BETA_Strat = case_when(ALT == A1 ~ BETA, REF == A1 ~ -1 * BETA))

  # Read in genotype matrix- this is the genotype matrix for ALT and mean center
  tp_file_path =  paste0(genos_prefix, i)
  X <- read_genos(geno_prefix = tp_file_path, betas = betas)

  # Compute PGS
  vecPGS[,i] <- pgs(X, betas)

  # Compute Va
  vecVa[i] <- calc_Va(X, betas)

}

# Combine chromosome
prs <- rowSums(matPGS)
Va <- sum(vecVa)


