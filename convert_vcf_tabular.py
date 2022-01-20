#! /usr/bin/env python
import sys, os, random 

if len(sys.argv) < 2: print("specify VCF file name",file=sys.stderr); sys.exit(); 

File  = open(sys.argv[1]);

###INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF|context|ancestral">

## A|stop_gained|HIGH|KMT2D|ENSG00000167548|Transcript|ENST00000301067|protein_coding|34/54||ENST00000301067.7:c.10183C>T|ENSP00000301067.7:p.Gln3395Ter|10183|10183|3395|Q/*|Caa/Taa||1||-1|SNV|1|HGNC|7133|YES||CCDS44873.1|ENSP00000301067|KMT2D_HUMAN|Q6PIA1_HUMAN&Q59FG6_HUMAN&F8VWW4_HUMAN|UPI0000EE84D6|||hmmpanther:PTHR22884&hmmpanther:PTHR22884:SF324|||||||||||||||||||POSITION:0.612916817142169&ANN_ORF:6351.44&MAX_ORF:6351.44|||HC|TGC|G

var_types = {'missense_variant':0,'frameshift_variant':0,'stop_gained':0,'splice_acceptor_variant':0,'splice_donor_variant':0,'synonymous_variant':0}
transcript = 'ENST00000301067'

print('#chrom','position','reference','alternate','annotation','impact','exon',sep='\t')
for line in File: 
    if line[0] == '#': continue
    var =line.strip().split('\t')
    if var[6] != 'PASS': continue 

    info = var[7].split(';')
    for pair in info: 
        if pair.split('=')[0] == 'CSQ': annotation = pair.split('=')[1].split(',')

    alist = []
    for anno in annotation: 
        if 'ENST00000301067' in anno: alist.append(anno)

    alleles = var[4].split(',')
    for i in range(len(alist)):
        A = alist[i].split('|')
        allele = A[0] 
        if allele in alleles and A[1] in var_types:
                print (var[0],var[1],var[3],allele,A[1],A[14],A[15],A[8],end='\n',sep='\t')
                var_types[A[1]] +=1
                #print(A)
File.close();

#for v in var_types: print(v,var_types[v])
#print("total number of LoF variants:",var_types['stop_gained'] + var_types['frameshift_variant'] + var_types['splice_acceptor_variant'] + var_types['splice_donor_variant'])
