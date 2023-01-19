## MED 263: Searching for disease mutations using DNA sequence data and bioinformatics 

We will be using Unix/Linux command line utilities such as 'grep', 'awk', 'sort', 'cut' and 'uniq' in this practical. If you are not familiar with these, you can find some basic information at this webpage: https://rsh249.github.io/bioinformatics/unix_shell.html 

## 1. Analysis of gene-specific variants using computational tools
Loss-of-function (LoF) mutations in genes are expected to have a strong impact on gene function. In the lecture, we learned that LoF mutations in the MLL2 (also known as KMT2D) cause Kabuki syndrome, a severe multi-system childhood disease. In this exercise, we will analyze coding variants in this gene using computational tools that annotate variants with their predicted functional impact.

(i) The ExAC database (currently part of gnomAD) contains variant calls from exome sequencing of 65,000 human individuals. First, we will use the 'tabix' tool to download the portion of the ExAC VCF file that contains all variants in the KMT2D gene. Tabix is a useful command line tool that works with tabular data (VCF files, bed files) to extract the subset of lines that overlap a genomic interval (start and end of the KMT2D gene in this example).

```Shell
tabix -h ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz 12:49412758-49453557  > KMT2D.ExAc.vcf
```

(ii) The VCF file stores information about the variant position, alleles and the impact of the variant on the coding or non-coding sequence
of the gene. We will use a python script to convert the VCF file into a simplified tabular format that contains the variant annotation relative to the canonical transcript for the KMT2D gene. 

```Shell
python3.6 convert_vcf_tabular.py  KMT2D.ExAc.vcf > KMT2D.ExAc.simplified.csv
```

The csv file can be viewed in a text editor or loaded as a spreadsheet. 

Using simple grep commands, we can count the number of different types of variants in this VCF file. LoF variants are of three types: stop_gain, splice_acceptor/splice_donor and frameshift.

```Shell
grep -E "stop_gained|splice_acceptor_variant|splice_donor_variant|frameshift" KMT2D.ExAc.simplified.csv | wc -l 
grep "missense" KMT2D.ExAc.simplified.csv | wc -l 
```

(iii) 6 of the 12 LoF variants in the variant file are splice-site disrupting variants (splice_acceptor or splice_donor). The deleterious impact of such variants can sometimes be 'mitigated' by usage of cryptic splice sites located near the canonical site. The deep-learning tool, SpliceAI [Jaganathan et al. Cell 2019](https://doi.org/10.1016/j.cell.2018.12.015), has high accuracy for predicting the effect of splicing altering variants. The splice variant, 12-49422743-T-C, is present in four individuals in the latest version of the gnomAD database. Therefore, it is highly unlikely that this variant is a LoF variant. We will use the spliceAI web-portal (https://spliceailookup.broadinstitute.org) to evaluate this variant further.

What does the prediction from SpliceAI tell us about the functional impact of this variant? Can we say with confidence that this variant is not a LoF variant?


(iv) Missense variants in the KMT2D gene have also been found to be pathogenic for Kabuki syndrome. However, a significant fraction of individuals in the normal population are also carriers of missense variants in KMT2D. Polyphen2 and SIFT are two of the most widely used tools for prioritizing missense mutations. Polyphen2 categorizes missense mutations as 'probably damaging' (D), 'possibly damaging (P), and benign (B) while SIFT classifies them as 'deleterious' (D) and 'neutral' (N). The file 'KMT2D.Kabuki.missense_predictions.csv' contains 57 missense variants (from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6488975/) identified in Kabuki syndrome individuals. Similarly, the file 'KMT2D.ExAc.missense_predictions.csv' has the ExAc missense variants. The Polyphen2 and SIFT predictions for each variant are included as separate columns in the two files.

To evaluate the ability of these two methods to distinguish between pathogenic missense variants (found in Kabuki syndrome patients) and benign missense variants (found in the ExAc database), we will count the number of missense variants predicted to be benign/deleterious in the two datasets: 

```Shell
python3.6 count_missense.py DATA/practical-1/KMT2D.Kabuki.missense_predictions.csv
python3.6 count_missense.py DATA/practical-1/KMT2D.ExAc.missense_predictions.csv
``` 

(v) We will use Fisher's exact test to determine if the computational classification of missense variants by SIFT can predict the pathogenicity of these variants. For this, we will tabulate the counts for SIFT in a 2x2 table ([[a,b],[c,d]]) and use the fisher_exact function from Python scipy.stats:

```
python3.6
>>> from scipy.stats import fisher_exact
>>> fisher_exact([[a,b],[c,d]])
```

The fisher_exact function returns two values: (i) the odds ratio and (ii) p-value (two-tailed probability of the odds-ratio being different than 1). Is the p-value significant? What is the interpretation of the odds ratio? 


## 2. Variant filtering in multiple related individuals (rare disease)
In the lecture, we talked about how DNA sequencing of related individuals can be used to find the genetic cause of rare diseases that affect individuals in
a family. This requires prioritizing variants based on a combination of (i) sharing by affected individuals, (ii) population allele frequency and (iii) impact on
gene function. In this exercise, we will use variants identified from exome sequencing of four individuals from a single family with a phenotype of early-onset
glaucoma (eye disease) to search for the disease causing variants. The variants were called using the GATK HaplotypeCaller and annotated using the Annovar tool.
The variants and genotypes have been summarized in a tabular format in the file "genotypes.coding.csv" (folder: DATA/practical-2)
Each variant has been annotated for its impact on genes and the allele frequency for each variant has also been obtained using data from the ExAc database. The allele frequency information is provided for two populations: European (EUR_AF, column 13) and South Asian (SAS_AF, column 14).
This file can be loaded into a spreadsheet as well. The four individuals correspond to:

* the mother (affected), label S1 and column 6 in the table
* child 1 (affected), label S2  and column 7 in the table
* child 2 (unaffected), label S3 and column 8 in the table
* child 3 (affected), label S4 and column 9 in the table

We want to search for variants that are shared by the three affected individuals (and not by the one unaffected individual), are rare in the population and also
affect the protein sequence. We will use a population allele frequency threshold of 0.1% to filter out common variants. Using awk, the following bash command can be used to search for potential disease-causing variants:

```Shell
cat DATA/practical-2/genotypes.coding.csv | awk '{FS="\t";} {if ($6 == "0/1" && $7 == "0/1" && $8 == "0/0" && $9 == "0/1" && $12 != "synonymous SNV" && $13 < 0.001 && $14 < 0.001) print; }' > candidates.csv
```

How many candidate variants are identified using the filter? 

We can use the Unix 'uniq' utility to find the number of distinct candidate genes with the variants:

```Shell
cut -f11 candidates.csv | uniq > candidates.genes
```

Next, we will determine if any of the candidate genes are already known to cause human diseases. For this, we will use the OMIM (https://omim.org) database that contains information about human genes and phenotypes. The file "omim.genes.disease" contains human genes and the corresponding phenotypes associated with them. We can search for each candidate gene in this table: 

```Shell
grep -w -f candidates.genes DATA/practical-2/omim.genes.disease 
```

Is there a gene that is known to cause an eye-related phenotype?


## 3. Using gene expression data for analyzing disease genes

Gene expression information can be used to prioritize genes for association with disease. The GTEx project (http://gtexportal.org/home/) has generated RNA-seq data
using more than 50 different tissues and cell-lines from hundreds of individuals. Summary data (RPKM values per gene for each tissue) is available for download from the GTEX website. We will use this data to analyze gene expression in disease-associated genes. 

File with RPKM values for all ENSEMBL transcripts and 50+ tissues/cell lines in a table format: 
* GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct (folder DATA/practical-2)

The first line of this file gives information about the tissues/cell-lines and each subsequent line has the expression information for an individual transcript. This file can easily be loaded into excel as well. 

(i) Extract the gene expression values for KMT2D from the data:

```Shell
grep -w KMT2D DATA/practical-2/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct
```
KMT2D is expressed at a high level across virtually all tissues which is consistent with the multi-organ phenotype associated with Kabuki syndrome. A visual plot of the RPKM values can be seen at http://gtexportal.org/home/gene/KMT2D or [here](DATA/practical-2/kmt2d_exp.png)


(ii) Compare the expression pattern for KMT2D to a gene RFX6 (discussed in Tuesday's lecture) which is expressed in a few tissues (stomach, pancreas, adrenal glands): http://gtexportal.org/home/gene/RFX6 or [here](DATA/practical-2/RFX6-expression.png)

```Shell
grep -w RFX6 DATA/practical-2/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct 
```

(iii) MLL2/KMT2D is the primary gene that is mutated in Kabuki syndrome (discussed in lecture). KDM6A is another gene that has been implicated in Kabuki syndrome (5-10% of individuals). This suggests that the genes should have a similar expression profile. We will calculate the correlation between the expression profiles of KMT2D and KDM6A using the scipy.stats.spearmanr function:

```Shell
python3.6 corr.py KMT2D KDM6A
```

The correlation between the expression values of the two genes is high. This is also apparent from visual inspection of plots for the two genes: http://gtexportal.org/home/gene/KMT2D and http://gtexportal.org/home/gene/KDM6A 



## Homework exercises

1. In class, we analyzed the predictions from the SIFT tool for missense variants. Perform the same analysis for the PolyPhen2 predictions. For this, you can group the 'probably damaging' (D) and 'possibly damaging (P) predictions into a single category. Which method (SIFT or Polyphen2) has a lower p-value? What can we infer about the predictive ability of these methods for classifying missense variants in the KMT2D gene? 


---

2. In the second practical exercise ("Variant filtering"), we used the population allele frequencies from ExAc to filter out common variants. Population allele frequencies for some variants can vary across populations, particularly for rare variants. The input file contains allele frequencies for each variant in two populations (European: column 13 and South Asian: column 14). The individuals sequenced in this study are from South Asia. Determine the number of candidate disease variants obtained by filtering using the allele frequency for each of the two population separately (use an allele frequency threshold of 0.001). Which population results in a lower number of candidate variants? Is this consistent with the ancestry of the individuals? Explain.

---

3. It is known that certain genes are expressed at very low levels in a specific cell type despite being expressed at a high level across almost all other cell types. Aberrant expression of genes in the "disallowed" cell type can result in disease phenotype(s). The file "disallowed.genes" contains a list of 20 genes that are known to be disallowed in pancreatic beta-cells. Use the GTEX RNA-seq expression data (used in class) to compare the expression of these genes in the tissue "pancreas" with the expression in other tissues. For each gene, calculate the ratio of the expression level in pancreas to the mean expression level in all other cell-types. For how many genes is this ratio less than 0.1? Use the OMIM database (https://www.omim.org) to determine if mutations in any of these genes are known to cause human disease that is related to the pancreas.







