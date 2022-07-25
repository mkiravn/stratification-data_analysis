## This script calculates lambda_T

## Requires: .eigenvec .eigenval .tvec

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_Tm.R <eigenvecs> <eigenvals> <tvec file> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

vecs_file = args[1] # eigenvectors
vals_file = args[2] # eigenvalues
tvec_file = args[3] # test vec file
out_file = args[4] # Name for out file

# Load test eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,2:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)

# Load test vector
TV <- fread(tvec_file)
Tvec_lat <- TV$latitude
Tvec_long <- TV$longitude

# Load Eigenvalues
vals <- fread(vals_file)
vals <- vals$V1

# Calculate Lambda T
lambda_T_lat <- t(Tvec_lat) %*% vecs %*% diag(vals) %*% t(vecs) %*% Tvec_lat
lambda_T_long <- t(Tvec_long) %*% vecs %*% diag(vals) %*% t(vecs) %*% Tvec_long

# Create output
out <- as.data.frame(cbind(lambda_T_lat, lambda_T_long))
colnames(out) <- c("latitude", "longitude")

# Save Lambda T
fwrite(out, out_file,row.names=F,quote=F,sep="\t", col.names = T)
