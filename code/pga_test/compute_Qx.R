## compute Qx
## This script reads in ascertained SNPs and the test vector and computes Qx and an empirical p-value

args=commandArgs(TRUE)

if(length(args)<9){stop("Rscript compute Qx.R <ascertained uncorrected> <ascertained lat> <ascertained long> <test panel prefix> <test vec file> <lambdaT file> <qx outfile> <pgs outfile>   ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

snp_u = args[1]
snp_lat = args[2]
snp_long = args[3]
snp_pc =args[4]
genos_prefix = args[5]
tvec_file = args[6]
lambdaT_file = args[7]
outfile_qx = args[8]
outfile_pgs = args[9]

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


main <- function(infile) {

  # Read in betas
  betas <- fread(infile)
  allBetasIDs <- betas$ID
  allBetas <- betas$BETA_Strat

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
  for (i in 1:num) {
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

  print(out)
  return(out)
}

# Compute Results
out <- matrix(NA, nrow = 4, ncol =4)
out[1, ] <- main(snp_u)
out[2, ] <- main(snp_lat)
out[3, ] <- main(snp_long)
out[4, ] <- main(snp_pc)

# Save output
colnames(out) <- c("Qx-Lat", "Qx-Long", "P-Lat", "P-Long")
row.names(out) <- c("Uncorrected", "TGWAS-Lat", "TGWAS-Long")
print(out)
fwrite(out, outfile_qx,row.names=T,quote=F,sep="\t", col.names = T)


# Function to just output PGS
main2 <- function(infile) {

  # Read in betas
  betas <- fread(infile)
  allBetasIDs <- betas$ID
  allBetas <- betas$BETA_Strat

  # Read in Genotypes
  X <- read_genos(genos_prefix, allBetasIDs)

  # Compute PGS
  sscore <- pgs(X, allBetas)

  return(sscore)
}


# Output File with all the PGS
fam <- fread(paste0(genos_prefix, ".psam"))
fam <- fam[,1:2]
fam$uncorrected <- main2(snp_u)
fam$lat <- main2(snp_lat)
fam$long <- main2(snp_long)
fam$PC <- main2(snp_pc)
fam$Tvec_Lat <- TV$latitude
fam$Tvec_Long <- TV$longitude

# Save output
print(head(fam))
fwrite(fam, outfile_pgs,row.names=F,quote=F,sep="\t", col.names = T)



