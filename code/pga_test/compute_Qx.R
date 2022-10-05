## compute Qx
## This script reads in ascertained SNPs and the test vector and computes Qx and an empirical p-value

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript compute Qx.R <prefix for ascertained snps> <test panel prefix> <test vec file> <lambdaT file> <qx outfile> <pgs outfile>   ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

snp_prefix = args[1]
genos_prefix = args[2]
tvec_file = args[3]
lambdaT_file = args[4]
outfile_qx = args[5]
outfile_pgs = args[6]

chr_start=1
chr_end=22
num = 1000


# Function to read in genotype matrix for a set of variants (counted allele is the Alternate)
read_genos <- function(geno_prefix, betas_id) {

  print(geno_prefix)
  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas_id
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=T) # Check this
  colnames(X) <- var.ids

  return(X)
}

# Function to calculate Va
calc_Va <- function(geno_mat, es) {

  # Get allele frequency in test panel
  freq <- colMeans(geno_mat) /2

  # Pull out effect sizes only
  effect_size <- es

  # Compute Va
  Va <- 2 * sum((effect_size)^2 * freq * (1 - freq))
  return(Va)
}


# Function to compute PGS
pgs <- function(X, betas) {

  # Mean center genotypes
  X <- apply(X, 2, function(y) y - mean(y))

  # Comput PGS
  Bhat_strat <- betas
  Z_strat <- X %*% Bhat_strat

  # Format output
  out <- cbind(Z_strat)
  colnames(out) <- c("STRAT")

  return(out)
}

# Function to calculate Qx
calc_Qx <- function(pgs, tvec, Va, lambda_T) {

  # Compute Qx Strat
  Ztest <- t(tvec) %*% pgs
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va*lambda_T)

  return(Qx_strat)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(betas, tvec, Va, X, lambda_T) {

  # Flip effect sizes
  b <- flip(betas)

  # Calculate PGS
  prs <- pgs(X, b)

  # Calculate Qx
  Qx <- calc_Qx(prs, tvec, Va, lambda_T)

  return(Qx)
}



# Load in Test vectors
TV <- fread(tvec_file)

# Load in LambdaT
lambdaT <- fread(lambdaT_file)


main <- function(beta_suffix) {

  # Read in all estimated betas by looping through chromosomes
  allBetasIDs <- c()
  allBetas<- c()
  for (i in chr_start:chr_end) {

    # Read in betas
    betas <- fread(paste0(snp_prefix, i, beta_suffix))

    # Deal with no betas by makes zeros
    if (nrow(betas) != 0) {

       # Flip betas to get the effect size of the ALT allele
       betas <- betas %>% mutate(BETA_Strat = case_when(ALT == A1 ~ BETA, REF == A1 ~ -1 * BETA))

       # Add betas to list
       allBetas <- c(allBetas, betas$BETA_Strat)
       allBetasIDs <- c(allBetasIDs, betas$ID)
   }
  }

  # Deal with no SNPs
  if (sum(allBetas) == 0) {
     out <- c(NA, NA, NA, NA)
  }

  else {
       # Read in Genotypes
       X <- read_genos(genos_prefix, allBetasIDs)

       # Compute PGS
       sscore <- pgs(X, allBetas)

       # Compute Va
       Va <- calc_Va(X, allBetas)

       # Compute Qx
       Qx_lat <- calc_Qx(sscore, TV$latitude,  Va, lambdaT$latitude)
       Qx_long <- calc_Qx(sscore, TV$longitude,  Va, lambdaT$longitude)

       # Generate Empirical null - Lat
       redraws <- matrix(0, ncol = 1, nrow = num)
       for (i in 1:num){
       	   redraws[i,] <- en(allBetas, TV$latitude, Va, X, lambdaT$latitude)
       }

       # Calculate empirical p-values
       all_strat <- redraws[,1]
       p_strat_en_lat <- length(all_strat[all_strat > Qx_lat[1,1]])/length(all_strat)

       # Generate Empirical null - Long
       redraws <- matrix(0, ncol = 1, nrow = num)
       for (i in 1:num){
       	   redraws[i,] <- en(allBetas, TV$longitude, Va, X, lambdaT$longitude)
       }

       # Calculate empirical p-values
       all_strat <- redraws[,1]
       p_strat_en_long <- length(all_strat[all_strat > Qx_long[1,1]])/length(all_strat)

       # Combine results
       out <- c(Qx_lat[1,1], Qx_long[1,1], p_strat_en_lat, p_strat_en_long)
       }

  print(out)
  return(out)
}

# Compute Results
out <- matrix(NA, nrow = 3, ncol =4)
out[1, ] <- main("_v3.Height.betas")
out[2, ] <- main("_v3.Height-Lat.betas")
out[3, ] <- main("_v3.Height-Long.betas")

# Save output
colnames(out) <- c("Qx-Lat", "Qx-Long", "P-Lat", "P-Long")
row.names(out) <- c("Uncorrected", "TGWAS-Lat", "TGWAS-Long")
print(out)
fwrite(out, outfile_qx,row.names=T,quote=F,sep="\t", col.names = T)


# Function to just output PGS
main2 <- function(beta_suffix) {

  # Read in all estimated betas by looping through chromosomes
  allBetasIDs <- c()
  allBetas<- c()
  for (i in chr_start:chr_end) {

    # Read in betas
    betas <- fread(paste0(snp_prefix, i, beta_suffix))

    # Deal with no betas by makes zeros
    if (nrow(betas) != 0) {

       # Flip betas to get the effect size of the ALT allele
       betas <- betas %>% mutate(BETA_Strat = case_when(ALT == A1 ~ BETA, REF == A1 ~ -1 * BETA))

       # Add betas to list
       allBetas <- c(allBetas, betas$BETA_Strat)
       allBetasIDs <- c(allBetasIDs, betas$ID)
   }
  }

  # Deal with no SNPs
  if (sum(allBetas) == 0) {
     sscore <- rep(NA, nrow(fam))
  }

  else {
       # Read in Genotypes
       X <- read_genos(genos_prefix, allBetasIDs)

       # Compute PGS
       sscore <- pgs(X, allBetas)
  }
  return(sscore)
}


# Output File with all the PGS
fam <- fread(paste0(genos_prefix, ".psam"))
fam <- fam[,1:2]
fam$uncorrected <- main2("_v3.Height.betas")
fam$lat <- main2("_v3.Height-Lat.betas")
fam$long <- main2("_v3.Height-Long.betas")
fam$Tvec_Lat <- TV$latitude
fam$Tvec_Long <- TV$longitude

# Save output
print(head(fam))
fwrite(fam, outfile_pgs,row.names=F,quote=F,sep="\t", col.names = T)



