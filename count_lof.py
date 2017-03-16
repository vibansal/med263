#! /usr/bin/env python
import sys, os, random 
from scipy import stats

## specify gene-name, get chromosome interval, get ExAc file, count loss-of-function mutations...


if len(sys.argv) < 2: print >>sys.stderr,"specify VCF file name"; sys.exit(); 

File  = open(sys.argv[1]);

combined_AC = 0;

for line in File: 
	if line[0] == '#': continue;
	var =line.strip().split('\t');
	if var[6] != 'PASS': continue; 
	alleles = var[4].split(',');
	AC1 = var[7].split(';')[0][3:]; AC = AC1.split(',');
	for i in xrange(len(alleles)):
		a = alleles[i]; 
		delta = len(a) -len(var[3]); 
		if ',' + a +'|stop_gained' in var[7]: 
			combined_AC += int(AC[i]); 
			print var[0],var[1],var[3],a,AC[i],'stop_gained';
		elif delta != 0 and '|frameshift' in var[7]: 
			print var[0],var[1],var[3],a,AC[i],'frameshift'
			combined_AC += int(AC[i]); 
		elif ',' + a +'|splice_donor' in var[7]: 
			combined_AC += int(AC[i]); 
			print var[0],var[1],var[3],a,AC[i],'splice_donor';
		elif ',' + a +'|splice_acceptor' in var[7]: 
			combined_AC += int(AC[i]); 
			print var[0],var[1],var[3],a,AC[i],'splice_acceptor';

File.close();

print combined_AC;
