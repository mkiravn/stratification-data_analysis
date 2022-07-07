# Snakefile to run UKBB  data analysis
CHR =["22"]
ROOT = ["/gpfs/data/berg-lab/jgblanc/stratification-data_analysis"]

rule all:
    input:
        expand("{root}/data/hgdp/variant_freq/hgdp_wgs.20190516.full.chr{chr}.afreq", root=ROOT,  chr=CHR)

## UKBB Genotype data processing

rule UKBB_begen_to_plink2:
    input:
        bgen="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr{chr}_v3.bgen",
        sample="/gpfs/data/berg-lab/data/ukbb/ukb27386_imp_v3_s487324.sample",
	pheno_ID="{root}/data/phenotypes/StandingHeight_50_IDs.txt"
    output:
        psam="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.psam",
	pvar="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pvar",
	ppgen="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen"
    params:
        prefix="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3"
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
	--threads 8 \
	--memory 38000 \
	--out {params.prefix} 
	"""

rule UKBB_freq:
    input:
        psam="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.psam",
        pvar="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pvar",
        ppgen="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen"
    output:
        freq="{root}/data/ukbb/variant_freq/ukb_imp_chr{chr}_v3.afreq"
    params:
        prefix_out="{root}/data/ukbb/variant_freq/ukb_imp_chr{chr}_v3",
	prefix_in="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3"
    shell:
        """
        plink2 --pfile {params.prefix_in} --freq \
        --threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """

## HGDP genotype data processing

rule HGDP_make_plink2:
    input:
        psam="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.pvar",
        ppgen="/gpfs/data/berg-lab/data/HGDP/plink2-files-hg19/hgdp_wgs.20190516.full.chr{chr}.pgen"
    output:
        psam="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pvar",
        ppgen="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pgen"
    params:
        prefix_out="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}",
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

rule HGDP_freq:
    input:
        psam="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pvar",
        ppgen="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pgen"
    output:
        freq="{root}/data/hgdp/variant_freq/hgdp_wgs.20190516.full.chr{chr}.afreq"
    params:
        prefix_in="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}",
        prefix_out="{root}/data/hgdp/variant_freq/hgdp_wgs.20190516.full.chr{chr}"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --freq \
	--threads 8 \
        --memory 38000 \
        --out {params.prefix_out}
        """



## Data processing for EUR analysis 

# Get separate list of mainland IDs and test sample IDs (FINs for now)
rule subset_1KG_IDs_EUR:
    input:
        samples="{root}/1000G_Genotypes/data/igsr_samples.tsv"
    output:
        target="{root}/1000G_Genotypes/FIN_IDs.txt",
	mainland="{root}/1000G_Genotypes/Mainland_IDs.txt"
    shell:
        """
	awk '$4 == "FIN"' {input.samples} | cut -f 1 | awk '{{print 0,$0}}' > {output.target}
	awk '$4 == "CEU" || $4 == "GBR" || $4 == "IBS" || $4 == "TSI"' {input.samples} | cut -f 1 | awk '{{print 0,$0}}' > {output.mainland}
	"""

# Get a list of SNPs with  
rule list_TP_SNPIDs:
    input:
        fam= "{root}/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim",
	bim="{root}/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim",
	bed="{root}/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bed"
    output:
        SNP_IDs="{root}/1000G_Genotypes/TestPanel_SNPIDs.txt",
        fam= "{root}/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR_maf0.05.fam",
        bim="{root}/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR_maf0.05.bim",
        bed="{root}/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR_maf0.05.bed"	
    params:
        prefix = "/scratch/jgblanc/stratification-data_analysis/1000G_Genotypes/data/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR"
    shell:
        """ 
	plink --bfile {params.prefix} --maf 0.05 --make-bed --out {params.prefix}_maf0.05
        cut -f2 {params.prefix}_maf0.05.bim > {output.SNP_IDs}
        """        


