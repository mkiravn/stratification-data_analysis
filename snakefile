# Snakefile to run UKBB  data analysis
CHR = []
for i in range(21, 23):
  CHR.append(str(i))
print(CHR)

rule all:
    input:
        "/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_genos.psam"

rule convert_bgen_to_plink:
    input:
        bgen="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr{chr}_v3.bgen",
        sample="/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb17346_imp_chr17_v3_s487378.sample"
    output:
        psam="/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{chr}_v3.psam",
        pvar="/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{chr}_v3.pvar",
        pgen="/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{chr}_v3.pgen"
    shell:
        """
	plink2 --bgen {input.bgen} ref-first --sample {input.sample} --maf 0.05 --rm-dup exclude-all --make-pgen --out /scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{wildcards.chr}_v3
	echo "{output.psam}" | cut -d'.' -f1 >> /scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/chromosome_list.txt
	"""

rule merge_chromosomes:
    input:
        psam=expand("/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{chr}_v3.psam", chr=CHR),
        pvar=expand("/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{chr}_v3.pvar", chr=CHR), 
	pgen=expand("/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_chr{chr}_v3.pgen", chr=CHR),
	chr_list="/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/chromosome_list.txt" 
    output:
        "/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_genos.psam",
        "/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_genos.pvar",
        "/scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_genos.pgen"
    shell:
        """
	plink2 --pmerge-list {input.chr_list} --out /scratch/jgblanc/stratification-data_analysis/UKBB_plink-files/ukb_imp_genos
	"""


