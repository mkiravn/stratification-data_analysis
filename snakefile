# Snakefile to run UKBB  data analysis
CHR =["22"]
ROOT = ["/scratch/jgblanc/stratification-data_analysis"]

#rule all:
#    input:
#        "/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_genos.psam"

rule list_UKBB_SNPIDs_chromosome:
    input:
        bgen="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr{chr}_v3.bgen",
        sample="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb17346_imp_chr17_v3_s487378.sample"
    output:
        ID="{root}/UKBB/variant_IDs/SNPs_{chr}.txt"
    params:
        prefix="/scratch/jgblanc/stratification-data_analysis/UKBB/temp"
    shell:
        """
	plink2 --bgen {input.bgen} ref-first --sample {input.sample} --maf 0.05 --rm-dup exclude-all --make-bgen --out {params.prefix} 
	cut -f 2 {params.prefix}.bim > {output.ID}
	rm {params.prefix}.*
	"""


## Data processing for Sardinia analysis

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


