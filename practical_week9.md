# MED 263: Searching for disease mutations using DNA sequence data and bioinformatics 


## 1. Analysis of loss-of-function mutations
Loss-of-function (LoF) mutations in genes have the strongest impact on gene function. In the lecture on Tuesday, we learned that LoF mutations in the MLL2 (also known as KMT2D)
cause Kabuki syndrome, a severe pediatric disease. Therefore, LoF mutations in this gene should be depleted in normal individuals. In this exercise, we will
estimate the frequency of LoF mutations in KMT2D in the ExAc database using command line tools and python.

(i) We will use the 'tabix' tool to download the portion of the ExAc VCF file that contains mutations in the KMT2D gene. 'tabix' is a very useful command line tool that works with tabular data (VCF files, bed files) to extract the subset of lines that overlap a genomic interval.

```Shell
tabix -h ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz 12:49412758-49453557  > KMT2D.ExAc.vcf
```

(ii) Using simple grep commands, we can identify the number of LoF variant sites in this VCF file. LoF variants are of three types: stop_gain, splice_acceptor/splice_donor and frameshift.

```Shell
grep "stop_gain" KMT2D.ExAc.vcf | grep PASS | wc -l
grep "splice_acceptor" KMT2D.ExAc.vcf | grep PASS | wc -l
```

What is the total number of LoF variant sites in the KMT2D gene ? Does the number of LoF variants match up with the ExAc website: http://exac.broadinstitute.org/gene/ENSG00000167548 (select LoF box)

(iii) Notice that some of the LoF sites are multi-allelic, i.e. the same base has multiple variant alleles. This information is represented in the VCF file on a single line but makes it difficult to parse it. We will use the python script "count_lof.py" to calculate the combined frequency of LoF variants in this gene.

```Shell
python count_lof.py KMT2D.ExAc.vcf
```

(iv) We will use LoF constraint scores (pLI) from the ExAc exome database to prioritize LoF mutations in an individual. The list of LoF mutations in an individual's exome has already been extracted from the VCF file.

```Shell
cat fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt | sort -k 2,2 > 1
cat NA12878.exome.vcf.annotated.LOFgenes | sort -k 1,1 > 2
join -1 2 -2 1 1 2 | cut -d ' ' -f1,20 > genes.LoF.constraintscores
```

What is the most constrained gene (high pLI score) ? Does this gene have a disease association in humans ?

## Homework exercise: Calculate the number of missense mutations in Exon 51 of KMT2D (12:49416372-49416658) in the ExAc database. Calculate the per-base rate of missense mutations in Exon 51 using the length of the exon. Do the same analysis for Exon 48 (12:49419964-49421105). Is the rate of missense mutations in Exon 51 higher than Exon 48 ?



## 2. Prioritizing disease genes using gene expression data

Gene expression information can be used to prioritize genes for association with disease independent of mutations. The GTEx project has generated RNA-seq data on more than 20 different tissues and cell-lines on > 50 individuals. Summary data (RPKM values per gene for each tissue) is available from the GTEX website (see file 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct'). We will use this data to analyze gene expression in the context of Kabuki disease.

(i) Extract the gene expression values for KMT2D from the data. KMT2D is expressed at a high level across virtually all tissues which is consistent with the multi-organ phenotype associated with Kabuki syndrome.

```Shell
grep KMT2D DATA/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct
```

(ii) Compare this expression pattern to a tissue-specific gene such as SLC30A8 which is expressed at a high level in the pancreas (RPKM = 3.9) and at very low levels in all other tissues.

```Shell
grep SLC30A8 DATA/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct
```

(iii) MLL2/KMT2D is the primary gene that is mutated in Kabuki syndrome (discussed in lecture). KDM6A is another gene that has been implicated in Kabuki syndrome. This suggests that the genes should have a similar expression profile. We will calculate the correlation between the expression profiles of KMT2D and KDM6A using the scipy.stats.spearmanr function:

```Shell
python corr.py KMT2D KDM6A
```

(iv) Correlation analysis can be used to find genes that have a similar expression profile to KMT2D.

```Shell
python corr.py KMT2D all > KMT2D.highcorrgenes
```

(v) Sort the list by the correlation coefficient value to find the top three genes whose expression is highly correlated with the expression profile of KMT2D. Use the constrained LoF scores data to determine if these genes are also constrained against LoF mutations.

```Shell
sort -k 3,3g KMT2D.highcorrgenes | tail -n 5
grep KMT2B fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt | cut -f2,20
grep BRPF1 fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt | cut -f2,20
```

### Homework exercise: Use the GTEX RNA-seq expression data to find genes that show a tissue-specific expression profile, i.e. genes for which the expression in the tissue with the maximum RPKM value is at least 5 times the RPKM values in all other tissues.



## 3. Phasing of heterozygous variants from sequence data.
In the lecture, we discussed how sequence reads can be used to infer haplotypes. In this exercise, we will use

Copy number variants result in decreased or increased read depth in aligned sequence data. Detecting copy number variants from exome sequencing is challenging due to the high variability in read depth across exons even in the absence of CNVs.




