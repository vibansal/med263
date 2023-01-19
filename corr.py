#! /usr/bin/env python
import sys, os, random 
from scipy import stats

## code to calculate correlation coefficient between gene expression levels of two genes (user specified) 

MIN_THRESH=0.6
if len(sys.argv) < 3: 
    print("specify two gene names",file=sys.stderr) 
    sys.exit() 

File = open('DATA/practical-2/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct');
if not File: print("gene expression file is missing"); sys.exit()
gene1 = sys.argv[1]; gene2 = sys.argv[2]; 

genenames = []; exps = [];

for line in File: 
	if line[0] == '#': continue;
	gene =line.strip().split();
	if len(gene) < 10: continue; 
	if gene[1] == gene1: exp1 = [float(i) for i in gene[2:]]; 
	elif gene[1] == gene2: exp2 = [float(i) for i in gene[2:]]; 
	else:  genenames.append(gene[1]); exps.append([float(i) for i in gene[2:]]); 


#print gene1,exp1; print gene2,exp2;

if gene2 != "all": print (gene1,gene2,stats.spearmanr(exp1,exp2))

else: 

	for i in range(len(genenames)):
		stat = stats.spearmanr(exp1,exps[i]); 
		if stat[0] >= MIN_THRESH: print(gene1,genenames[i],stat[0],stat[1])

