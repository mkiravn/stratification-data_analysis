## Overlapping SNPs
## This script takes two sets of .afreq files and the SDS file and nd subsets each to variants over 1% and then gets the intersection of the three lists

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript overlapping_snps.R <ukbb.freq> <test panel.freq> <sds pol file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

ukbb_file = args[1]
tp_file = args[2]
sds_file = args[3]
outfile = args[4]

# Load freqency files
ukbb <- fread(ukbb_file)
tp <- fread(tp_file)
sds <- fread(sds_file) # if this is the polarised file then we don't have to worry about flipping bc it already matches ukbb

# Subset to SNPs > 1% MAF
ukbb <- subset(ukbb, ukbb$ALT_FREQS > 0.01 & ukbb$ALT_FREQS < 0.99)
tp <- subset(tp, tp$ALT_FREQS > 0.01 & tp$ALT_FREQS < 0.99)
sds <- subset(sds, sds$DAF > 0.01 & sds$DAF < 0.99)
# Subset to SNPs with missingness rate < 5%
# no need to do this with SDS
tp <- subset(tp, tp$OBS_CT > (0.95 * max(tp$OBS_CT)))
ukbb <- subset(ukbb, ukbb$OBS_CT > (0.95 *  max(ukbb$OBS_CT)))

# Get overlapping SNPs with same alt/ref
sds <- sds %>% rename("#CHROM"="CHR")
matched <- inner_join(ukbb, tp, by = c("#CHROM", "ID", "ALT", "REF")) %>%
              inner_join(sds,by=c("#CHROM", "ID", "ALT", "REF"))
matched <- matched %>% select(`#CHROM`,ID,REF,ALT )

# Get overlapping SNPs with opposite alt/ref
tmp1 <-inner_join(ukbb,sds,by=c("#CHROM", "ID")) %>%
          inner_join(tp, by = c("#CHROM", "ID")) 
flipped <- subset(tmp1, tmp1$REF.x == tmp1$ALT.y & tmp1$ALT.x == tmp1$REF.y)
flipped <- flipped %>% select("#CHROM", ID, REF.x, REF.y)
colnames(flipped) <- c("#CHROM", "ID", "REF", "ALT")

# Get complete list of SNPs and list the reference allele in UKBB
out <- rbind(matched, flipped) %>% select(ID, REF)
fwrite(out,outfile, row.names = F, col.names = T, quote = F, sep = "\t")

