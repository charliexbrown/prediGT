#!/usr/bin/bash
set -euo pipefail

##define flags and args of the script
while getopts p:e:o: flag
do
    case "${flag}" in
        p) path=${OPTARG};; ##path where to find input files
        e) ext=${OPTARG};;  #extension of input files
        o) out=${OPTARG};; ##output file name 
    esac
done

echo $path
echo $ext
echo $out

##create out_dir if not existing 
if [ -d "$path/results" ]; then
    :
else
    mkdir -p "$path/results"
fi

out_dir="$path/results"

## extract significant gene in each tissue predixcan output & list gene - tissue

for i in "$path"/*"$ext"
do
   filename=$(basename "$i")
   tex_name=$(echo "$filename" | sed 's/^[^_]*_\([^\.]*\)\.csv$/\1/')
   #tex_name=$(echo $i|sed 's/.csv//')
   echo "$tex_name"
   awk -F "," -v OFS="\t" -v tex="$tex_name" '{if($5<0.05) print $2, tex}' "$i" >> ${out_dir}/sig_genes_tex.csv ##extract p<0.05 genes 
   awk -F "," -v OFS="\t" '{if ($5<0.05) print $1, $2}' "$i" >> ${out_dir}/sig_ens_gene.csv 
done 

sort ${out_dir}/sig_ens_gene.csv | uniq > ${out_dir}/sig_ens_gene_dic.csv

##groupby per gene & list tissues
awk '{arr[$1] = arr[$1] "," $2} END {for (i in arr) print i "\t" substr(arr[i], 2)}' ${out_dir}/sig_genes_tex.csv | sort -k1 > ${out_dir}/gene_tex_list.csv

## count of tissue per gene
cut -f 1 ${out_dir}/sig_genes_tex.csv | sort | uniq -c |  sed 's/^[ \t]*//' | sed 's/ /\t/' > ${out_dir}/sig_genes_counts.csv

## merge count file & tissue list
join -1 1 -2 2 ${out_dir}/gene_tex_list.csv ${out_dir}/sig_genes_counts.csv | sed 's/ /\t/g' | sort -k3 -n -r > ${out_dir}/${out}

#mv ${out}  "${out_dir}/" 
#mv sig_ens_gene_dic.csv "${out_dir}/"

#rm sig_genes_counts.csv
#rm sig_genes_tex.csv
#rm gene_tex_list.csv
#rm sig_ens_gene.csv
