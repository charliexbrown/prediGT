# prediGT

Pipeline for Genome-Wide Association Study (GWAS) and Transcriptome-Wide Association Study (TWAS) analyses. 

## HOW IT WORKS
prediGT is divided into two modules. 
First it tests association between genotypes with (a) trait(s) of interest using plink. Input genotyping data should be provided as plink bfiles for each chromosome (*chr.bed* + *chr.bim* + *chr.fam*). Phenotypes file should be also be provided following plink format (i.e. FID, IID, phenotype(s)).
Second, it investigates the association between genetically regulated gene expression (i.e. eQTLs) and the trait, using MetaXcan tools (https://github.com/hakyimlab/MetaXcan). Pre-computed predictive models are downloaded from GTEx v8 project to run TWAS.

## INSTALLATION
1. Clone this repository
2. Create conda environment using *snakemake.yaml*
```
 conda env create -f snakemake.yaml
 conda activate snakemake 
```

3. install MetaXcan software (https://github.com/hakyimlab/MetaXcan) 
```
 git clone https://github.com/hakyimlab/MetaXcan
```
4. run script down_gtex.sh to download predictive models from GTEx (https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/) and polish file names.

## HOW TO RUN
**Run from snakemake directory.**
1. Edit *config.yaml* file with paths and parameters 

2. run:
```
snakemake -s Snakemake --configfile config.yaml 
```

Notes:
Add ``` --latency-wait n ``` if pipeline struggles to create output files.


