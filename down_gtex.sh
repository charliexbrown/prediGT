#!/bin/sh

mkdir gtex_models
cd gtex_models

## download files
wget https://zenodo.org/record/3518299/files/mashr_eqtl.tar
tar -xvf mashr_eqtl.tar

wget https://zenodo.org/record/3518299/files/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz
cd  eqtl/mashr

for i in mashr*; do new_name=$(echo $i | sed 's/mashr_//g'); mv $i $new_name; done
