## Format covar
## This script reads in all the potential covariates (age at reccruitment, genetic sex, array type, TGWAS) and formats them for plink

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript conccat_Tm.R <prefix to Tm chromosomes> <chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

fam = fread(args[1])
TGWAS = fread(args[2])
aar = fread(args[3])
sex = fread(args[4])
array = fread(args[5])
outfile = args[6]

# Format
df <- as.data.frame(matrix(NA, nrow = nrow(TGWAS), ncol = 2))
colnames(df) <- c("#FID", "IID")
df$`#FID` <- fam$`#FID`
df$IID <- fam$IID

# Merge Age
colnames(aar) <- c("V1", "IID", "Age", "V2")
m1 <- inner_join(df, aar, by = c("IID" = "IID")) %>% select("#FID", "IID", "Age")

# Merge sex
colnames(sex) <- c("V1", "IID", "Sex", "V2")
m2 <- inner_join(m1, sex, by = c("IID" = "IID")) %>% select("#FID", "IID", "Age", "Sex")

# Merge array type
colnames(array) <- c("V1", "IID", "Array", "V2")
array <- array %>% mutate(Array =  case_when(Array < 0 ~ 0, Array > 0 ~ 1))
m3 <- inner_join(m2, array, by = c("IID" = "IID")) %>% select("#FID", "IID", "Age", "Sex", "Array")

# Merge TGWAS
out <- cbind(m3, TGWAS$latitude, TGWAS$longitude)
colnames(out) <- c("#FID", "IID", "Age", "Sex", "Array", "latitude", "longitude")

# Drop missing values
out <- drop_na(out)

# Write to output file
fwrite(out, outfile, row.names = F, col.names = T, quote = F, sep = "\t")



