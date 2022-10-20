# Snakefile to run HGDP (full data) / UKBB  data analysis
CHR =[]
for i in range(1, 23):
  CHR.append(str(i))
ROOT = ["/gpfs/data/berg-lab/jgblanc/stratification-data_analysis"]
DATASET = ["EUR"]
PVAL = ["p_1"]


def get_params(x):
  out = x.split("_")[1]
  return out

def get_size_minus_one(x):
  if  x == "ALL":
    out = 928
  elif x == "EUR":
    out = 154
  else:
    out = "ERROR"
  return out

rule all:
    input:
        expand("{root}/data/ukbb-hgdp/pga_test/{dataset}/{pval}/Qx.txt", root=ROOT, chr=CHR, dataset = DATASET, pval=PVAL)


## UKBB Genotype data processing

rule UKBB_begen_to_plink2:
    input:
        bgen="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr{chr}_v3.bgen",
        sample="/gpfs/data/berg-lab/data/ukbb/ukb22828_c22_b0_v3_s487192.sample", ## Need to fix sample file
	pheno_ID="data/phenotypes/StandingHeight_50_IDs.txt"
    output:
        psam="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.psam",
	pvar="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pvar",
	pgen="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen"
    params:
        prefix="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3"
    shell:
        """
	plink2 --bgen {input.bgen} ref-first \
	--sample {input.sample} \
	--keep {input.pheno_ID} \
	--maf 0.01 \
	--rm-dup exclude-all \
	--snps-only \
	--max-alleles 2 \
	--make-pgen \
	--set-all-var-ids @:# \
	--threads 20 \
	--memory 38000 \
	--out {params.prefix}
	"""

rule UKBB_freq:
    input:
        psam="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.psam",
        pvar="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pvar",
        pgen="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen"
    output:
        freq="{root}/data/ukbb/variant_freq/ukb_imp_chr{chr}_v3.afreq"
    params:
        prefix_out="{root}/data/ukbb/variant_freq/ukb_imp_chr{chr}_v3",
	      prefix_in="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3"
    shell:
        """
        plink2 --pfile {params.prefix_in} --freq \
        --threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """

## HGDP genotype data processing

rule get_EUR_sample_list:
    input:
        "{root}/data/hgdp/hgdp_wgs.20190516.metadata.txt"
    output:
        "{root}/data/hgdp/hgdp_wgs.20190516.EUR_samples.txt"
    shell:
        """
        Rscript code/genotypes/get_EUR_samples.R {input} {output}
        """

rule HGDP_make_plink2_all:
    input:
        psam="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.pgen"
    output:
        psam="{root}/data/hgdp/plink2-files/ALL/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/ALL/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/hgdp/plink2-files/ALL/hgdp_wgs.20190516.full.chr{chr}.pgen"
    params:
        prefix_out="{root}/data/hgdp/plink2-files/ALL/hgdp_wgs.20190516.full.chr{chr}",
	      prefix_in="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --maf 0.01 \
        --rm-dup exclude-all \
        --snps-only \
        --max-alleles 2 \
        --make-pgen \
        --set-all-var-ids @:# \
        --threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """

rule HGDP_make_plink2_EUR:
    input:
        psam="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.pgen",
        samples="{root}/data/hgdp/hgdp_wgs.20190516.EUR_samples.txt"
    output:
        psam="{root}/data/hgdp/plink2-files/EUR/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/EUR/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/hgdp/plink2-files/EUR/hgdp_wgs.20190516.full.chr{chr}.pgen"
    params:
        prefix_out="{root}/data/hgdp/plink2-files/EUR/hgdp_wgs.20190516.full.chr{chr}",
	      prefix_in="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --maf 0.01 \
        --keep {input.samples} \
        --rm-dup exclude-all \
        --snps-only \
        --max-alleles 2 \
        --make-pgen \
        --set-all-var-ids @:# \
        --threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """

rule HGDP_freq:
    input:
        psam="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pgen"
    output:
        freq="{root}/data/hgdp/variant_freq/{dataset}/hgdp_wgs.20190516.full.chr{chr}.afreq"
    params:
        prefix_in="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}",
        prefix_out="{root}/data/hgdp/variant_freq/{dataset}/hgdp_wgs.20190516.full.chr{chr}"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --freq \
	--threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """

## Get overlapping set of final SNPS

# Right now set for 1%
rule get_overlapping_snps:
    input:
        freq_hgdp="{root}/data/hgdp/variant_freq/{dataset}/hgdp_wgs.20190516.full.chr{chr}.afreq",
	freq_ukbb="{root}/data/ukbb/variant_freq/ukb_imp_chr{chr}_v3.afreq"
    output:
        "{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    shell:
        """
        Rscript code/genotypes/overlapping_snps.R {input.freq_ukbb} {input.freq_hgdp} {output}
        """

## Recode HGDP with UKBB ref/alt alleles and save to new directory

rule HGDP_recode:
    input:
        psam="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pgen",
      	snp_list="{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    output:
        psam="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pgen"
    params:
        prefix_in="{root}/data/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}",
        prefix_out="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --extract {input.snp_list} \
	--ref-allele {input.snp_list} \
	--make-pgen \
        --out {params.prefix_out}
        """

## Make latitude and longitude test vector

rule make_Tvec_cordinates:
    input:
        psam="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr22.psam",
        populations="{root}/data/hgdp/hgdp_wgs.20190516.metadata.txt"
    output:
        "{root}/data/ukbb-hgdp/calculate_Tm/ALL/Tvec_cordinates.txt",
        "{root}/data/ukbb-hgdp/calculate_Tm/EUR/Tvec_cordinates.txt"
    shell:
        """
        Rscript code/calculate_Tm/make_Tvec_hgdp_cordinates.R {input.psam} {input.populations} {wildcards.root}/data/ukbb-hgdp/calculate_Tm/ALL/Tvec {wildcards.root}/data/ukbb-hgdp/calculate_Tm/EUR/Tvec
        """

## Compute TGWAS

# Each chromosome individually
rule project_Tvec_chr:
    input:
        Tvec="{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/Tvec_cordinates.txt",
        tp_genos="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}.pgen",
        gp_genos="/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen",
        overlap_snps="{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    output:
        "{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/Tm_{chr}.txt"
    params:
        tp_prefix = "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.chr{chr}",
        gp_prefix = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3",
        tvec_prefix = "{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/Tvec",
        out_prefix = "{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/"
    shell:
        """
        Rscript code/calculate_Tm/project_Tvec_chr.R {params.tp_prefix} {params.gp_prefix} {params.tvec_prefix} {params.out_prefix} {input.overlap_snps} {output}
        """

# Add together individal chromosomes
rule concat_chr_Tm:
    input:
        expand("{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/Tm_{chr}.txt", chr = CHR, root=ROOT, dataset = DATASET)
    output:
        "{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/TGWAS.txt"
    params:
        chromosomes = CHR
    shell:
        """
        Rscript code/calculate_Tm/concat_Tm.R {wildcards.root}/data/ukbb-hgdp/calculate_Tm/{wildcards.dataset}/Tm {output} {params.chromosomes}
        """

## Run GWAS

rule format_covars:
    input:
        fam = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr22_v3.psam",
        TGWAS = "{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/TGWAS.txt",
        aar = "/gpfs/data/berg-lab/data/ukbb/phenotypes/age_at_recruitment_21022.txt",
        sex = "/gpfs/data/berg-lab/data/ukbb/phenotypes/genetic_sex_22001.txt",
        array = "/gpfs/data/berg-lab/data/ukbb/phenotypes/genotype_measurement_batch_22000.txt",
        PCs = "/gpfs/data/berg-lab/data/ukbb/phenotypes/genetic_PC_22009.txt"
    output:
        "{root}/data/ukbb-hgdp/run_gwas/{dataset}/covars.txt"
    shell:
        """
        Rscript code/run_gwas/format_covars.R {input.fam} {input.TGWAS} {input.aar} {input.sex} {input.array} {input.PCs} {output}
        """

rule run_gwas_uncorrected:
    input:
        genos = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen",
        covar = "{root}/data/ukbb-hgdp/run_gwas/{dataset}/covars.txt",
        pheno = "{root}/data/phenotypes/StandingHeight_50.txt",
	snp_list = "{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    output:
        "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3.Height.glm.linear"
    params:
        pfile = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3",
        out = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3"
    shell:
        """
        plink2 \
        --pfile {params.pfile} \
	--extract {input.snp_list} \
        --glm hide-covar \
        --covar {input.covar} \
        --covar-col-nums 3-5 \
        --pheno {input.pheno} \
        --pheno-name Height \
        --out {params.out}
        """

rule run_gwas_latitude:
    input:
        genos = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen",
        covar = "{root}/data/ukbb-hgdp/run_gwas/{dataset}/covars.txt",
        pheno = "{root}/data/phenotypes/StandingHeight_50.txt",
	snp_list = "{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    output:
        "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-Lat.Height.glm.linear"
    params:
        pfile = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3",
        out = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-Lat"
    shell:
        """
        plink2 \
        --pfile {params.pfile} \
	--extract {input.snp_list} \
        --glm hide-covar \
        --covar {input.covar} \
        --covar-col-nums 3-6 \
        --pheno {input.pheno} \
        --pheno-name Height \
        --out {params.out}
        """

rule run_gwas_longitude:
    input:
        genos = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen",
        covar = "{root}/data/ukbb-hgdp/run_gwas/{dataset}/covars.txt",
        pheno = "{root}/data/phenotypes/StandingHeight_50.txt",
	      snp_list = "{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    output:
        "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-Long.Height.glm.linear"
    params:
        pfile = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3",
        out = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-Long"
    shell:
        """
        plink2 \
        --pfile {params.pfile} \
        --extract {input.snp_list} \
        --glm hide-covar \
        --covar {input.covar} \
        --covar-col-nums 3-5,7 \
        --pheno {input.pheno} \
        --pheno-name Height \
        --out {params.out}
        """

rule run_gwas_PCs:
    input:
        genos = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen",
        covar = "{root}/data/ukbb-hgdp/run_gwas/{dataset}/covars.txt",
        pheno = "{root}/data/phenotypes/StandingHeight_50.txt",
	      snp_list = "{root}/data/ukbb-hgdp/variants/{dataset}/snps_chr{chr}.txt"
    output:
        "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-PCs.Height.glm.linear"
    params:
        pfile = "/scratch/jgblanc/ukbb/plink2-files/ukb_imp_chr{chr}_v3",
        out = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-PCs"
    shell:
        """
        plink2 \
        --pfile {params.pfile} \
        --extract {input.snp_list} \
        --glm hide-covar \
        --covar {input.covar} \
        --covar-col-nums 3-5,8-47 \
        --pheno {input.pheno} \
        --pheno-name Height \
        --out {params.out}
        """


## Ascertain SNPs for PGS

rule ascertain_snps:
    input:
        betas_uncorrected = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3.Height.glm.linear",
        betas_lat = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-Lat.Height.glm.linear",
        betas_long= "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-Long.Height.glm.linear",
        betas_PCs = "{root}/data/ukbb-hgdp/run_gwas/effect_sizes/{dataset}/ukb_imp_chr{chr}_v3-PCs.Height.glm.linear",
        block = "{root}/data/LD_blocks/fourier_ls-all_parsed.bed"
    output:
        snps_uncorrected = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_chr{chr}_v3.Height.betas",
        snps_lat = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_chr{chr}_v3.Height-Lat.betas",
        snps_long = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_chr{chr}_v3.Height-Long.betas",
        snps_PCs = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_chr{chr}_v3.Height-PCs.betas"
    params:
        pt = lambda wildcards: get_params(wildcards.pval)
    shell:
        """
        Rscript code/run_gwas/pick_snps.R {input.block} {input.betas_uncorrected} {input.betas_lat} {input.betas_long} {input.betas_PCs} {output.snps_uncorrected} {output.snps_lat} {output.snps_long} {output.snps_PCs} {params.pt}
        """

## Do PGA Test

rule combine_hgdp_chromosomes:
    input:
        hgdp_genos = expand("{{root}}/data/ukbb-hgdp/hgdp/plink2-files/{{dataset}}/hgdp_wgs.20190516.full.chr{chr}.psam", chr = CHR)
    output:
        "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.psam",
        "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.pvar",
        "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.pgen"
    params:
        out_prefix = "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full",
        file_list = expand("{{root}}/data/ukbb-hgdp/hgdp/plink2-files/{{dataset}}/hgdp_wgs.20190516.full.chr{chr}", chr = CHR)
    shell:
        """
        echo {params.file_list} > tmp1.txt
        cat tmp1.txt | tr ' ' '\n' > chrs.txt
        plink2 \
        --pmerge-list chrs.txt \
        --out {params.out_prefix}
        rm tmp1.txt
        rm chrs.txt
        """

rule PCA:
    input:
        psam="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.psam",
        pvar="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.pvar",
        pgen="{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.pgen"
    output:
        "{root}/data/ukbb-hgdp/hgdp/PCA/{dataset}/pca.eigenvec",
        "{root}/data/ukbb-hgdp/hgdp/PCA/{dataset}/pca.eigenval"
    params:
        out_prefix = "{root}/data/ukbb-hgdp/hgdp/PCA/{dataset}/pca",
        pfile_prefix = "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full",
        n_minus_1 = lambda wildcards: get_size_minus_one(wildcards.dataset)
    shell:
        """
        plink2 \
        --pfile {params.pfile_prefix} \
        --pca {params.n_minus_1} \
        --out {params.out_prefix}
        """

rule calc_lambdaT:
    input:
        vecs="{root}/data/ukbb-hgdp/hgdp/PCA/{dataset}/pca.eigenvec",
        vals="{root}/data/ukbb-hgdp/hgdp/PCA/{dataset}/pca.eigenval",
        tvec="{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/Tvec_cordinates.txt",
    output:
        "{root}/data/ukbb-hgdp/pga_test/{dataset}/Lambda_T.txt"
    shell:
        """
      	Rscript code/pga_test/calc_lambdaT.R {input.vecs} {input.vals} {input.tvec} {output}
	      """

rule concat_ascertained_snps:
    input:
        snps_uncorrected = expand("{{root}}/data/ukbb-hgdp/run_gwas/ascertained/{{dataset}}/{{pval}}/ukb_imp_chr{chr}_v3.Height.betas", chr = CHR),
        snps_lat = expand("{{root}}/data/ukbb-hgdp/run_gwas/ascertained/{{dataset}}/{{pval}}/ukb_imp_chr{chr}_v3.Height-Lat.betas", chr=CHR),
        snps_long = expand("{{root}}/data/ukbb-hgdp/run_gwas/ascertained/{{dataset}}/{{pval}}/ukb_imp_chr{chr}_v3.Height-Long.betas", chr = CHR),
        snps_PCs = expand("{{root}}/data/ukbb-hgdp/run_gwas/ascertained/{{dataset}}/{{pval}}/ukb_imp_chr{chr}_v3.Height-PCs.betas", chr = CHR)
    output:
        uncorrected = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height.betas",
        lat = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height-Lat.betas",
        long = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height-Long.betas",
        PCs = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height-PCs.betas"
    params:
        snp_prefix = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_chr"
    shell:
        """
        Rscript code/run_gwas/concat_snps.R {params.snp_prefix} {output.uncorrected} {output.lat} {output.long} {output.PCs}
        """


rule compute_Qx:
    input:
        snps_uncorrected = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height.betas",
        snps_lat = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height-Lat.betas",
        snps_long = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height-Long.betas",
        snps_PCs = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_all_v3.Height-PCs.betas",
        Tvec = "{root}/data/ukbb-hgdp/calculate_Tm/{dataset}/Tvec_cordinates.txt",
        hgdp_genos = "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full.psam",
        lambdaT = "{root}/data/ukbb-hgdp/pga_test/{dataset}/Lambda_T.txt"
    output:
        Qx = "{root}/data/ukbb-hgdp/pga_test/{dataset}/{pval}/Qx.txt",
        PGS = "{root}/data/ukbb-hgdp/pga_test/{dataset}/{pval}/PGS.txt"
    params:
        snp_prefix = "{root}/data/ukbb-hgdp/run_gwas/ascertained/{dataset}/{pval}/ukb_imp_chr",
        tp_prefix = "{root}/data/ukbb-hgdp/hgdp/plink2-files/{dataset}/hgdp_wgs.20190516.full"
    shell:
        """
        Rscript code/pga_test/compute_Qx.R {input.snps_uncorrected} {input.snps_lat} {input.snps_long} {input.snps_PCs} {params.tp_prefix} {input.Tvec} {input.lambdaT} {output.Qx} {output.PGS}
        """

