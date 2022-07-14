#!/bin/bash
pfile_path=$1
beta_path=$2
outfile=$3

plink2 \
  --pfile $pfile_path \
  --score $beta_path header-read center cols=dosagesum,scoresums \
  --out $outfile
