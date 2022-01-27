## MED 263: Searching for disease mutations using DNA sequence data and bioinformatics 

We will be using Unix/Linux command line utilities such as 'grep', 'awk', 'sort', 'cut' and 'uniq' in this practical. If you are not familiar with these, you can find some basic information at this webpage: https://rsh249.github.io/bioinformatics/unix_shell.html 

## 1. Prioritizing disease genes using population variant data
Loss-of-function (LoF) mutations in genes are expected to have a strong impact on gene function. In the lecture, we learned that LoF mutations in the MLL2 (also known as KMT2D) cause Kabuki syndrome, a severe multi-system childhood disease. Therefore, LoF mutations in this gene should be depleted in normal individuals. In this exercise, we will determine the frequency of LoF mutations in KMT2D in the ExAc database (65,000 individuals) using command line tools and python.

(i) First, we will use the 'tabix' tool to download the portion of the ExAc VCF file that contains all mutations in the KMT2D gene. 'tabix' is a very useful command line tool that works with tabular data (VCF files, bed files) to extract the subset of lines that overlap a genomic interval (start and end of the KMT2D gene in this example).

```Shell
tabix -h ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz 12:49412758-49453557  > KMT2D.ExAc.vcf
```

(ii) The VCF file stores information about the variant position, alleles and the impact of the variant on the coding or non-coding sequence
of the gene. We will use a python script to convert the VCF file into a simplified tabular format that contains the variant annotation relative to the canonical transcript for the KMT2D gene. 

```Shell
python3.6 convert_vcf_tabular.py  KMT2D.ExAc.vcf > KMT2D.ExAc.simplified.csv
```

The csv file can be viewed in a text editor or loaded as a spreadsheet. 

Using simple grep commands, we can count the number of LoF variant sites in this VCF file. LoF variants are of three types: stop_gain, splice_acceptor/splice_donor and frameshift.

```Shell
grep -E "stop_gained|splice_acceptor_variant|splice_donor_variant|frameshift" KMT2D.ExAc.simplified.csv | wc -l 
```

What is the total number of LoF variant sites in the KMT2D gene in the ExAc database?  

(iii) The ExAc database provides "LoF" constraint scores (pLI score) for human genes based on the observed and expected counts of LoF mutations in each gene. The constraint scores range from 0 (no constraint) to 1 (completely constrained). For the KMT2D gene, the observed:expected ratio is 10:126 which implies that the observed number of LoF variants is much less than expected. We will use the list of scores to find the rank of the KMT2D (MLL2) gene. The data file "fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt" contains the summary of the constraint scores. 

```Shell
sort -k 2,2 DATA/practical-1/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt | cut -f 2,20 > allgenes.constraint.scores
cat allgenes.constraint.scores | sort -k 2,2gr | awk '{ a += 1; if ($1 == "KMT2D") print $1,a; }'
```
You can also load this file into excel and sort to find the rank.


(iv) Missense variants in the KMT2D gene have also been identified in individuals with Kabuki syndrome. We will use the simplified variant file to count the number of missense and silent (synonymous) variants in the KMT2D gene: 

```Shell
grep -E "missense_variant" KMT2D.ExAc.simplified.csv | wc -l 
grep -E "synonymous_variant" KMT2D.ExAc.simplified.csv | wc -l 
```

What is the ratio of the number of missense variants and the number of synonymous variants? Exon 38 has been observed to be enriched in disease-causing mutations in this gene. Therefore, we would expect to see a lower than expected number of missense variants in the normal population in this exon relative to synonymous variants. To check this, let us calculate the number of missense and synonymous variants in Exon 38: 


```Shell
grep -E "missense_variant" KMT2D.ExAc.simplified.csv | awk '{ if ($8 ==  "38/54") print; }' | wc -l 
grep -E "synonymous_variant" KMT2D.ExAc.simplified.csv | awk '{ if ($8 ==  "38/54") print; }' | wc -l 
```

Is the missense_variant:synonymous_variant ratio for Exon 38 smaller than that for the entire gene? Is the difference statistically significant?

This type of analysis can be used to prioritize specific regions of genes that are depleted of missense variants in normal individuals and are likely to harbor disease-causing mutations. For more details, you can read this paper: https://www.biorxiv.org/content/10.1101/148353v1


## 2. Prioritizing disease genes using gene expression data 

Gene expression information can be used to prioritize genes for association with disease. The GTEx project (http://gtexportal.org/home/) has generated RNA-seq data
using more than 50 different tissues and cell-lines from hundreds of individuals. Summary data (RPKM values per gene for each tissue) is available for download from the GTEX website. We will use this data to analyze gene expression in disease-associated genes. 

File with RPKM values for all ENSEMBL transcripts and 50+ tissues/cell lines in a table format: 
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct (folder DATA/practical-2)

The first line of this file gives information about the tissues/cell-lines and each subsequent line has the expression information for an individual transcript. This file can easily be loaded into excel as well. 

(i) Extract the gene expression values for KMT2D from the data:

```Shell
grep KMT2D DATA/practical-2/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct
```
KMT2D is expressed at a high level across virtually all tissues which is consistent with the multi-organ phenotype associated with Kabuki syndrome. A visual plot of the RPKM values can be seen at http://gtexportal.org/home/gene/KMT2D or [here](DATA/practical-2/kmt2d_exp.png)


(ii) Compare the expression pattern for KMT2D to a gene RFX6 (discussed in Tuesday's lecture) which is expressed in a few tissues (stomach, pancreas, adrenal glands): http://gtexportal.org/home/gene/RFX6 or [here](DATA/practical-2/RFX6-expression.png)

```Shell
grep RFX6 DATA/practical-2/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct 
```

(iii) MLL2/KMT2D is the primary gene that is mutated in Kabuki syndrome (discussed in lecture). KDM6A is another gene that has been implicated in Kabuki syndrome. This suggests that the genes should have a similar expression profile. We will calculate the correlation between the expression profiles of KMT2D and KDM6A using the scipy.stats.spearmanr function:

```Shell
python3.6 corr.py KMT2D KDM6A
```

The correlation between the expression values of the two genes is high. This is also apparent from visual inspection of plots for the two genes: http://gtexportal.org/home/gene/KMT2D and http://gtexportal.org/home/gene/KDM6A 


(iv) Next, we will use correlation analysis to find genes that have a very similar (corr. coefficient > 0.9) expression profile to KMT2D.

```Shell
python3.6 corr.py KMT2D all > KMT2D.highcorrgenes
```

(v) We can sort the list of genes by the correlation coefficient value to find the top three genes whose expression is highly correlated with the expression profile of KMT2D. Using the constrained LoF scores data, we will determine if these genes are also constrained against LoF mutations.

```Shell
sort -k 3,3gr KMT2D.highcorrgenes | head
grep KMT2B allgenes.constraint.scores
grep BRPF1 allgenes.constraint.scores
```

These genes correspond to epigenetic regulators or histone-modifying proteins and have been linked to rare childhood diseases: [KMT2B](https://www.omim.org/entry/617284) and [BRPF1](https://www.omim.org/entry/617333)


## 3. Variant filtering in rare disease
In the lecture, we talked about how DNA sequencing of related individuals can be used to find the genetic cause of rare diseases that affect individuals in
a family. This requires prioritizing variants based on a combination of (i) sharing by affected individuals, (ii) population allele frequency and (iii) impact on
gene function. In this exercise, we will use variants identified from exome sequencing of four individuals from a single family with a phenotype of early-onset
glaucoma (eye disease) to search for the disease causing variants. The variants were called using the GATK HaplotypeCaller and annotated using the Annovar tool.
The variants and genotypes have been summarized in a tabular format in the file "genotypes.coding.csv". 
Each variant has been annotated for its impact on genes and the allele frequency for each variant has also been obtained using data from the ExAc database. The allele frequency information is provided for two populations: European (column 13) and South Asian (column 14).
This file can be loaded into a spreadsheet as well. The four individuals correspond to:

* the mother (affected), label S1 and column 6 in the table
* child 1 (affected), label S2  and column 7 in the table
* child 2 (unaffected), label S3 and column 8 in the table
* child 3 (affected), label S4 and column 9 in the table

We want to search for variants that are shared by the three affected individuals (and not by the one unaffected individual), are rare in the population and also
affect the protein sequence. We will use a population allele frequency threshold of 0.1% to filter out common variants. Using awk, the following bash command can be used to search for potential disease-causing variants:

```Shell
cat DATA/practical-3/genotypes.coding.csv | awk '{FS="\t";} {if ($6 == "0/1" && $7 == "0/1" && $8 == "0/0" && $9 == "0/1" && $12 != "synonymous SNV" && $13 < 0.001 && $14 < 0.001) print; }' > candidates.csv
```

How many candidate variants are identified using the filter? 

We can use the Unix 'uniq' utility to find the number of distinct candidate genes with the variants:

```Shell
cut -f11 candidates.csv | uniq > candidates.genes
```

Next, we will determine if any of the candidate genes are already known to cause human diseases. For this, we will use the OMIM (https://omim.org) database that contains information about human genes and phenotypes. The file "omim.genes.disease" contains human genes and the corresponding phenotypes associated with them. We can search for each candidate gene in this table: 

```Shell
grep -w -f candidates.genes DATA/practical-3/omim.genes.disease 
```

Is there a gene that is known to cause an eye-related phenotype?


## Homework exercises

1. Genes that are expressed primarily in a single tissue or cell-type are likely to be important for the function of that tissue and relevant for diseases that affect that specific tissue. 
Use the GTEX RNA-seq expression data to find genes that show a tissue-specific expression profile, i.e. genes for which the expression in the tissue with the maximum RPKM value is at least 5 times the RPKM values in all other tissues. Report the top 5 genes that are primarily expressed in 'pancreas'. 

---

2. It is known that certain genes are expressed at very low levels in a specific cell type despite being expressed universally across almost all other cell types. One example of such a cell type is pancreatic beta-cells that release insulin in response to glucose. The file "disallowed.genes" has a list of 20 genes that are known to be disallowed in pancreatic beta-cells. Use the GTEX RNA-seq expression data to compare the expression of these gene in pancreas with the expression in other tissues. For each gene, calculate the ratio of the expression level in pancreas to the mean expression level in all other cell-types. For how many genes is this ratio less than 0.1? Use the OMIM database to determine if mutations in any of these genes are known to cause a insulin secretion defect due to aberrant expression. 

---

3. For the third practical exercise ("Variant filtering"), we used the population allele frequencies from ExAc to filter out common variants. The input file contains allele frequencies for each variant in two population (European: column 13 and South Asian: column 14). Determine the number of candidate disease variants obtained by filtering using the allele frequency for each of the two population separately (threshold of 0.001). Which population results in a lower number of candidate variants? What can we infer about the ancestry of the sequenced individuals from this? 







