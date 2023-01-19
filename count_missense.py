#! /usr/bin/env python
import sys, os, random 

File  = open(sys.argv[1]);

transcript = 'ENST00000301067'

SIFT = {}
Polyphen = {}
n=0
for line in File: 
    if line[0] == '#': continue
    info = line.strip().split('\t')
    try: SIFT[info[5]] +=1
    except KeyError: SIFT[info[5]] =1
    try: Polyphen[info[6]] +=1
    except KeyError: Polyphen[info[6]] =1
    n +=1

File.close();

print('SIFT-counts',end=' ')
for s in sorted(SIFT.keys()): print (s,SIFT[s],round(SIFT[s]/n,3),end=' ',sep=':')
print('\n')


print('Polyphen-counts',end=' ')
for s in sorted(Polyphen.keys()): print (s,Polyphen[s],round(Polyphen[s]/n,3),end=' ',sep=':')
print()
