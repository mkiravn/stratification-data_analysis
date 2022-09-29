#!/bin/bash
pfile_path=$1
beta_path=$2
outfile=$3
overlap_snps=$4

plink2 \
  --pfile $pfile_path \
  --extract $overlap_snps \
  --score $beta_path header-read variance-standardize cols=dosagesum,scoresums \
  --out $outfile
