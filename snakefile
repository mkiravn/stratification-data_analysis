# Snakefile to run UKBB  data analysis
CHR =["21", "22"]
#for i in range(1, 23):
#  CHR.append(str(i))
ROOT = ["/gpfs/data/berg-lab/jgblanc/stratification-data_analysis"]

rule all:
    input:
        expand("{root}/data/ukbb-hgdp/calculate_Tm/Tm_{chr}.txt", root=ROOT,  chr=CHR)

## UKBB Genotype data processing

rule UKBB_begen_to_plink2:
    input:
        bgen="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr{chr}_v3.bgen",
        sample="/gpfs/data/berg-lab/data/ukbb/ukb27386_imp_v3_s487324.sample", ## Need to fix sample file
	pheno_ID="{root}/data/phenotypes/StandingHeight_50_IDs.txt"
    output:
        psam="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.psam",
	pvar="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pvar",
	pgen="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen"
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
        pgen="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen"
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
        pgen="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pgen"
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

## Get overlapping SNPS

rule get_overlapping_snps:
    input:
        freq_hgdp="{root}/data/hgdp/variant_freq/hgdp_wgs.20190516.full.chr{chr}.afreq",
	freq_ukbb="{root}/data/ukbb/variant_freq/ukb_imp_chr{chr}_v3.afreq"
    output:
        "{root}/data/ukbb-hgdp/variants/snps_chr{chr}.txt"
    shell:
        """
        Rscript code/genotypes/overlapping_snps.R {input.freq_ukbb} {input.freq_hgdp} {output}
        """

## Recode HGDP with UKBB ref/alt alleles

rule HGDP_recode:
    input:
        psam="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pgen",
	snp_list="{root}/data/ukbb-hgdp/variants/snps_chr{chr}.txt"
    output:
        psam="{root}/data/ukbb-hgdp/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.psam",
        pvar="{root}/data/ukbb-hgdp/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pvar",
        pgen="{root}/data/ukbb-hgdp/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pgen"
    params:
        prefix_in="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}",
        prefix_out="{root}/data/ukbb-hgdp/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}"
    shell:
        """
        plink2 --pfile {params.prefix_in} \
        --extract {input.snp_list} \
	--ref-allele {input.snp_list} \
	--make-pgen \
        --out {params.prefix_out}
        """

## Make latitdue test vector

rule make_Tvec_cordinates:
    input:
        psam="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr22.psam",
        populations="{root}/data/hgdp/hgdp_wgs.20190516.metadata.txt"
    output:
        "{root}/data/ukbb-hgdp/calculate_Tm/Tvec_cordinates.txt"
    shell:
        """
        Rscript code/calculate_Tm/make_Tvec_hgdp_cordinates.R {input.psam} {input.populations} {wildcards.root}/data/ukbb-hgdp/calculate_Tm/Tvec
        """

## Compute TGWAS

# Each chromosome individually
rule project_Tvec_chr:
    input:
        Tvec="{root}/data/ukbb-hgdp/calculate_Tm/Tvec_cordinates.txt",
        tp_genos="{root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{chr}.pgen",
        gp_genos="{root}/data/ukbb/plink2-files/ukb_imp_chr{chr}_v3.pgen",
        overlap_snps="{root}/data/ukbb-hgdp/variants/snps_chr{chr}.txt"
    output:
        "{root}/data/ukbb-hgdp/calculate_Tm/Tm_{chr}.txt"
    shell:
        """
        Rscript code/calculate_Tm/project_Tvec_chr.R {wildcards.root}/data/hgdp/plink2-files/hgdp_wgs.20190516.full.chr{wildcards.chr} {wildcards.root}/data/ukbb/plink2-files/ukb_imp_chr{wildcards.chr}_v3 {wildcards.root}/data/ukbb-hgdp/calculate_Tm/Tvec {wildcards.root}/data/ukbb-hgdp/calculate_Tm/ {input.overlap_snps} {output}
        """

# Add together individal chromosomes
rule concat_chr_Tm:
    input:
        expand("{root}/data/ukbb-hgdp/calculate_Tm/Tm_{chr}.txt", chr = CHR)
    output:
        "{root}/data/ukbb-hgdp/calculate_Tm/TGWAS.txt"
    params:
        chromosomes = CHR
    shell:
        """
        Rscript code/calculate_Tm/concat_Tm.R {wildcards.root}/data/ukbb-hgdp/calculate_Tm/Tm {params.chromosomes} {output}
        """






