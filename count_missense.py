#! /usr/bin/env python
import sys, os, random 

File  = open(sys.argv[1]);

transcript = 'ENST00000301067'

SIFT = {}
Polyphen = {}
for line in File: 
    if line[0] == '#': continue
    info = line.strip().split('\t')
    try: SIFT[info[5]] +=1
    except KeyError: SIFT[info[5]] =1
    try: Polyphen[info[6]] +=1
    except KeyError: Polyphen[info[6]] =1

File.close();

print('SIFT-counts',end=' ')
for s in SIFT.keys(): print (s,SIFT[s],end=' ',sep=':')
print()
print('Polyphen-counts',end=' ')
for s in Polyphen.keys(): print (s,Polyphen[s],end=' ',sep=':')
print()
