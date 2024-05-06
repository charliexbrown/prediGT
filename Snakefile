configfile: "config.yaml"
import glob
import pandas as pd

# Snakefile for GWAS and transcriptomic imputation analysis

# Define input files
bfile_path=config['bfile_path']
w_dir=config['working_dir']
metaxcan_dir=config['metaxcan_dir']
pheno_file=config['pheno_file']
covar_file=config['covar_file']
sig_gwas=config['sig_pval']
glm_type=config['glm_type']

# define wildcards
## Create list of chromosomes
chromosomes=[]
for i in os.listdir(bfile_path):
        if ".fam" in i:
                chromosomes.append(i.split(".")[0])

print("chr: ", chromosomes)

## Create  list of tissues
snake_folder= os.getcwd()
gtex_models=[]
model_file=os.path.join(w_dir, "model_db_list.txt")
model_path=os.path.join(snake_folder, "gtex_models", "eqtl", "mashr")
print(model_path)

with open(model_file, "w") as f:
        for i in os.listdir(model_path):
                if ".db" in i:
                        f.write("{}/{}\n".format(model_path, i))
                        gtex_models.append(i.split(".")[0])
print("\ntex:", gtex_models)

## Create a list of phenotypes from header of pheno.txt, to use as params in rule gwas
header=[]
with open(pheno_file, "r") as f:
        header= f.readline().strip().split("\t")
phenotype=header[2:]
print("\npheno:", phenotype)

# define output: gwas results, gwas plots, single-tissue predixcan results, combined single-tissue predixcan, multi-tissue multixcan
rule all:
    input:
        expand("{bfile_path}/{chromosomes}.{ext}", chromosomes=chromosomes, bfile_path=bfile_path, ext=['fam', 'bed', 'bim']),
        expand("{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas_results.txt", w_dir=w_dir, phenotype=phenotype),
        expand("{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas.txt", w_dir=w_dir, phenotype=phenotype),
        expand("{w_dir}/{phenotype}/gwas_results/{chromosomes}.{phenotype}.{glm_type}", w_dir=w_dir, phenotype=phenotype, chromosomes=chromosomes, glm_type=glm_type),
        expand("{w_dir}/{phenotype}/spredixcan_output/{phenotype}_{gtex_models}.csv", w_dir=w_dir, phenotype=phenotype, gtex_models=gtex_models),
        expand("{w_dir}/{phenotype}/spredixcan_output/results/multixcan/{phenotype}_multix.csv", w_dir=w_dir, phenotype=phenotype),
        expand("{w_dir}/{phenotype}/spredixcan_output/results/{phenotype}_spredixcan_alltex.csv", w_dir=w_dir, phenotype=phenotype),
        expand("{w_dir}/{phenotype}/plots/manhattan.png", w_dir=w_dir, phenotype=phenotype),
        expand("{w_dir}/{phenotype}/plots/qqplot.png", w_dir=w_dir, phenotype=phenotype)


##############
#### GWAS ####
##############
## function to define output files from either logistic or linear regression in plink2

def get_glm_type(wildcards, temp):
    if os.path.isfile(temp.format(chromosomes=wildcards.chromosomes, phenotype=wildcards.phenotype, glm_type="glm.linear")):
        return "glm.linear"
    elif os.path.isfile(temp.format(chromosomes=wildcards.chromosomes, phenotype=wildcards.phenotype, glm_type="glm.logistic")):
        return "glm.logistic"

# Define rule for GWAS analysis for each phenotype
rule gwas:
    input:
        expand("{bfile_path}/{chromosomes}.{ext}",  chromosomes=chromosomes, ext=['fam', 'bed', 'bim'], bfile_path=bfile_path),
        pheno_inp = pheno_file,
        covariate = covar_file if covar_file else []
    output:
        "{w_dir}/{phenotype}/gwas_results/{chromosomes}.{phenotype}.{glm_type}"
    params:
        temp="{w_dir}/{chromosomes}.{phenotype}.{glm_type}",
        out_dir="{w_dir}/{phenotype}/gwas_results/"
    run:
        glm_type = get_glm_type(wildcards, params.temp)
        if input.covariate:
            shell("plink2 --bfile {bfile_path}/{wildcards.chromosomes} --pheno {input.pheno_inp} --covar {input.covariate} --glm hide-covar --allow-extra-chr --out {w_dir}/{wildcards.chromosomes} --threads 6")
        else:
            shell("plink2 --bfile {bfile_path}/{wildcards.chromosomes} --pheno {input.pheno_inp} --glm --allow-extra-chr --out {w_dir}/{wildcards.chromosomes} --threads 6")
        shell("mv {params.temp} {params.out_dir}")

# Define rule to concatenate GWAS results for all chromosomes for each phenotype
rule concat_gwas:
    input:
        expand("{w_dir}/{phenotype}/gwas_results/{chromosomes}.{phenotype}.{glm_type}", w_dir=w_dir, phenotype=phenotype, chromosomes=chromosomes, glm_type=glm_type)
    output:
        "{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas_results.txt"
    shell:
        """
        head -1 {input[0]} > {output}
        for file in {input}; do tail -n +2 $file  >> {output}; done
        """

# Fix SNP_ID if not present (i.e. ./NA) and harmonize with SNP_ID in gtex models (chr_pos_alt_ref_b38)
rule fix_id:
    input:
        "{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas_results.txt"
    output:
        "{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas.txt"
    run:
        df = pd.read_csv(input[0], sep="\t")
        if "chr" in df.iloc[:,0].values:
            df["ID"] = df[df.columns[0]].astype(str) + '_' + df[df.columns[1]].astype(str) + '_' + df["ALT"] + '_' + df["REF"] + "_b38"
        else:
            df["ID"]="chr"+ df[df.columns[0]].astype(str) + '_' + df[df.columns[1]].astype(str) + '_' + df["ALT"] + '_' + df["REF"] + "_b38"
        if "BETA" in df.columns:
            pass
        else:
            df.rename(columns={"LOG(OR)_SE":"BETA"}, inplace=True)
        df.to_csv(output[0], sep="\t", index=None)

# Plot gwas results
rule plot_gwas:
    input:
        "{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas.txt"
    output:
        "{w_dir}/{phenotype}/plots/manhattan.png",
        "{w_dir}/{phenotype}/plots/qqplot.png"
    params:
        path="{w_dir}/{phenotype}/plots/",
        rscript=f"{snake_folder}/Scripts/gwas_plot.R"
    shell:
        """
	Rscript {params.rscript} {input} {params.path}
        """


# Extract significant snps in gwas. Used later when retrieving significant variants associated with significant association from spredixcan
rule sig_snps:
    input:
        "{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas.txt"
    output:
        "{w_dir}/{phenotype}/{phenotype}_sig_snps.csv"
    params:
        sig_gwas=sig_gwas
    shell:
       "awk -F '\t' '$13<{params.sig_gwas} {{print $1}}' {input} > {output}"


##########################################
#### TWAS: GENE EXPRESSION IMPUTATION ####
##########################################

# Define rule for transcriptomic imputation with SPrediXcan and SMultiXcan
rule spredixcan:
    input:
        genotypes="{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas.txt",
        model_db= f"{model_path}/{{gtex_models}}.db",
        covariance=f"{model_path}/{{gtex_models}}.txt.gz"
    output:
        "{w_dir}/{phenotype}/spredixcan_output/{phenotype}_{gtex_models}.csv"
    params:
        metaxcan_dir=metaxcan_dir,
    run:
        try:
            shell("python {params.metaxcan_dir}/SPrediXcan.py --gwas_file {input.genotypes} --model_db_path {input.model_db} --covariance {input.covariance} --snp_column ID --model_db_snp_key varID --effect_allele_column ALT --non_effect_allele_column REF --beta_column BETA --pvalue_column P --keep_non_rsid --additional_output --throw --output_file {output}")
        except: ## overcomes error in generating output file in case no association are found
            shell("python {params.metaxcan_dir}/SPrediXcan.py --gwas_file {input.genotypes} --model_db_path {input.model_db} --covariance {input.covariance} --snp_column ID --model_db_snp_key varID --effect_allele_column ALT --non_effect_allele_column REF --beta_column BETA --pvalue_column P --keep_non_rsid --throw --output_file {output}")

rule smultixcan:
    input:
        genotypes="{w_dir}/{phenotype}/gwas_results/{phenotype}_gwas.txt",
        model_db=expand("{model_path}/{gtex_models}.db", model_path=model_path, gtex_models=gtex_models),
        covariance="gtex_models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz",
        spredix=lambda wildcards: expand("{w_dir}/{phenotype}/spredixcan_output/{phenotype}_{gtex_models}.csv", w_dir=w_dir, phenotype=wildcards.phenotype, gtex_models=gtex_models)
    output:
        "{w_dir}/{phenotype}/spredixcan_output/results/multixcan/{phenotype}_multix.csv",
    params:
        meta_folder= lambda wildcards: "{w_dir}/{phenotype}/spredixcan_output/".format(w_dir=wildcards.w_dir, phenotype=wildcards.phenotype),
        db_path= model_path,
        metaxcan_dir=metaxcan_dir
    shell:
        "python {params.metaxcan_dir}/SMulTiXcan.py --models_folder {params.db_path} --models_name_pattern '(.*).db' --model_db_snp_key varID --gwas_file {input.genotypes} --snp_covariance {input.covariance} --snp_column ID --keep_non_rsid --effect_allele_column ALT  --non_effect_allele_column REF --beta_column BETA --se_column SE --pvalue_column P --metaxcan_folder {params.meta_folder}  --metaxcan_file_name_parse_pattern '(.*)_(.*).csv' --output {output} --throw --cutoff_condition_number 30 --verbosity 4"

# Post TWAS processing
## extract significant genes in all tissues of predixcan output
rule sig_predixcan:
    input:
        expand("{w_dir}/{phenotype}/spredixcan_output/{phenotype}_{gtex_models}.csv", w_dir=w_dir, phenotype=phenotype, gtex_models=gtex_models)
    output:
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/{phenotype}_merged_predixcan.csv"),
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/sig_genes_counts.csv"),
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/sig_genes_tex.csv"),
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/gene_tex_list.csv"),
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/sig_ens_gene.csv"),
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/sig_ens_gene_dic.csv")
    params:
        in_path="{w_dir}/{phenotype}/spredixcan_output",
        out_name="{phenotype}_merged_predixcan.csv",
        script=f"{snake_folder}/Scripts/extract_sig_gene.sh"
    shell:
        """
        chmod +x {params.script}
        {params.script} -p {params.in_path} -e ".csv" -o {params.out_name}
        """

## extract significant variants for each significant gene in predixcan
rule extract_variants:
    input:
        sig_ensid_list= "{w_dir}/{phenotype}/spredixcan_output/results/sig_ens_gene_dic.csv",
        sig_vars= "{w_dir}/{phenotype}/{phenotype}_sig_snps.csv",
        model_db_list= "{w_dir}/model_db_list.txt"
    output:
        temporary("{w_dir}/{phenotype}/spredixcan_output/results/sig_gene_rsid.csv")
    params:
        out_path="{w_dir}/{phenotype}/spredixcan_output/results/",
        script=f"{snake_folder}/Scripts/rsid_fromdb.py"
    run:
        try:
           shell("python {params.script} {input.sig_ensid_list} {input.sig_vars} {input.model_db_list} --out_path {params.out_path}")
        except:
           print("no significant variants associated to significant TWAS associations")
           shell(">{output}")

rule polish_final:
    input:
        spred_counts="{w_dir}/{phenotype}/spredixcan_output/results/{phenotype}_merged_predixcan.csv",
        gene_rs="{w_dir}/{phenotype}/spredixcan_output/results/sig_gene_rsid.csv"
    output:
        final="{w_dir}/{phenotype}/spredixcan_output/results/{phenotype}_spredixcan_alltex.csv"
    run:
        import pandas as pd
        try:
             spred=pd.read_csv(input.spred_counts, header=None, sep='\t', names=['gene_name', 'tissues', 'n_tix'])
             rsids=pd.read_csv(input.gene_rs, sep='\t')
             merged=spred.merge(rsids, on='gene_name')
             print(merged.head())
             merged.to_csv(output.final, index=False, sep='\t', na_rep="NA")
        except:
             shell("> {output}")
rule rm_log:
    input:
        expand("{w_dir}/{chromosomes}.log", w_dir=w_dir, chromosomes=chromosomes)
    shell:
        "rm {input}"
