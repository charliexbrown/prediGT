import pandas as pd 
import sqlite3
from multiprocessing import Pool
import multiprocessing 
import sys 
import os 
import time
import argparse

##define script arguments, start_time and cores
parser=argparse.ArgumentParser()
parser.add_argument("sig_genes", help="file with list of sig genes ensid and gene name")
parser.add_argument("sig_vars", help="file with list of significant snps")
parser.add_argument("model_list", help='file with list of models')
parser.add_argument("--out_path", help="optional path for output file", default='nd')

args=parser.parse_args()

if len(sys.argv) < 4:
    print("invalid arguments: \nusage is: {} sig_ensid_list sig_vars_list model_list optional_output_path".format(sys.argv[0]))
    
 
st = time.time()
print(multiprocessing.cpu_count())

#define input files for ensid_list, database_paths and gwas_file (to extract ids for significant variants) 
try:
	gene_ens=pd.read_csv(sys.argv[1], sep="\t", header=None)
except: 
	gene_ens=pd.read_csv(sys.argv[1], sep=",", header=None)
print(gene_ens)

gene_ens.columns=['ensid', 'gene_name']
gene_list=gene_ens.ensid.tolist()
    
database_paths=[]
with open(sys.argv[3], "r") as f:
      database_paths=f.read().splitlines()

sig_vars=[]
with open(sys.argv[2], "r") as f:
	sig_vars=f.read().splitlines()

#define output file
if args.out_path != 'nd':
	output_file=os.path.join(args.out_path, 'sig_gene_rsid.csv')
else:
	output_file='sig_gene_rsid.csv'

#define funct to extract rsids and varids from sqlite3 db
def extract_rsid(database_path, gene):
	conn = sqlite3.connect(database_path)
	c = conn.cursor() #create a cursor obj to execute sql queries on the db
	rsid_list = []
	varid_list=[]
	c.execute("SELECT rsid,varid FROM weights WHERE gene = ?", (gene,)) ##execute sql query on the db to extract rsid and varid for gene of interest
	result = c.fetchall() #retrieve query and store in result variable
	if result: #check for result and store rsid to list and varid to other list
		for row in result:
			rsid_list.append(row[0])
			varid_list.append(row[1])
	conn.close()
	return rsid_list, varid_list #return two lists

#loop over list of gene to extract all rsids annotated in all databases

list_sig_rsid=[]
list_sig_varid=[]
list_all_rsid=[]
list_all_varid=[]

left_genes=len(gene_list)

print (f"processing {left_genes} genes...")

for gene in gene_list:
	left_genes-=1
	print(f"{left_genes} genes left to process", end='\r')
	
	#create empty sets for storing all rsid and varids per gene (uniq values)
	all_rsid=set() ##all snps mapped to gene
	all_varid=set() 
	sig_rsid=set() ##significant snps mapped to gene
	sig_varid=set()
	
	# apply extract_rsid for each database to extract rsid and varid for the gene from each database
	ensid_lists=[]
	for database_path in database_paths:
		ensid_list=extract_rsid(database_path, gene)
		ensid_lists.append(ensid_list)
	# ensid_lists is a tuple with two elements corrisponding to rsids and varids: assign each element to a separate list	
	for rs,var in ensid_lists:
		for r,v in zip(rs,var):
			all_rsid.add(r)
			all_varid.add(v)
			if v in sig_vars: #check if variants are in the gwas_file to store only ids for significant variants 
				sig_rsid.add(r)
				sig_varid.add(v)
	# store all rsid lists and all varid lists in a whole list (to create df)
	list_sig_rsid.append(sig_rsid)
	list_sig_varid.append(sig_varid)
	list_all_rsid.append(all_rsid)
	list_all_varid.append(all_varid)

df=pd.DataFrame({'gene':gene_list, 'sig_rsids': list_sig_rsid, 'sig_varids': list_sig_varid, 'all_rsids': list_all_rsid, 'all_varids':list_all_varid} )
new=df.merge(gene_ens, left_on='gene', right_on='ensid').drop(['gene', 'ensid'], axis=1).sort_values(by='gene_name')

##refine df to substitute empty set with "NA"
prova=new.applymap(lambda x: "NA" if isinstance(x, set) and len(x) == 0 else x)
prova[['gene_name', 'sig_rsids', 'sig_varids', 'all_rsids', 'all_varids']].to_csv(output_file, sep='\t', index=None)


##get end time and print exec time
et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
