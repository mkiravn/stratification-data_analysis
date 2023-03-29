## Re-polarises SDS scores (AA/DA) to match UKBB (REF/AlT)

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript repolarise_SDS.R <sds file> <ukbb file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

sds_file = args[1]
ukbb_file = args[2]
outfile = args[3]

####################
## Functions #######
####################
####################
joined <- SDS %>% 
  rename("#CHROM"="CHR") %>%
  inner_join(ukbb,by=c("#CHROM","POS","ID"))

joined <- joined %>%
  rowwise() %>%
  mutate(SDS=ifelse(AA==REF & DA==ALT,SDS,-SDS))

fwrite(joined, outfile, row.names = F, col.names = T, quote = F, sep = "\t")
