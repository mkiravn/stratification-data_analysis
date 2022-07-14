#!/bin/bash
pfile_path=$1
pheno_path=$2
test_type=$3
outfile=$4

plink2 \
  --pfile $pfile_path \
  --glm omit-ref allow-no-covars \
  --pheno $pheno_path \
  --pheno-name $test_type \
  --geno-counts \
  --out $outfile
