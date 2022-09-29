#!/bin/bash
bgen_path=$1
beta_path=$2
outfile=$3
overlap_snps=$4
sample=$5
pheno_ID=$6

plink2 \
    --bgen $bgen_path ref-first \
    --sample $sample \
    --keep $pheno_ID \
    --set-all-var-ids @:# \
    --extract $overlap_snps \
    --score $beta_path header-read variance-standardize cols=dosagesum,scoresums \
    --out $outfile
