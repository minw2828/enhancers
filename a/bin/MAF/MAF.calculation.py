#!/usr/bin/env python

################################################################################
#
# Description:
# This module calculates minor allele frequency of the input genotype file.
#
# Documentation:
# http://biopython.org/DIST/docs/api/Bio.SearchIO.BlastIO-module.html
#
# Usage:
# python <MODULE_NAME> <INPUT_FILE> <OUTPUT_FILE> <INFMT>
#
# Author:
# Min Wang (min.wang@depi.vic.gov.au)
#
# Date Created:
# 09 April 2016
# 
# The program is designed to be used with Python 2.6 and 2.7.
#
# Date modified and reason:
#
################################################################################


import sys
sys.path.insert(0, '/group/dairy/Min/geno2pheno/software/python/Python-3.5.1/')
import resource
import os
from ggplot import *

INFILE1    = sys.argv[1] # genotype file, no rowname nor column names, sum row
INFILE2    = sys.argv[2] # map file 
INFILE3    = sys.argv[3] # enhancer lable file 
OUTFILE    = sys.argv[4]


def read_file(infile):
    f = open(infile)
    fc = [line.strip() for line in f.readlines()]
    f.close()
    return(fc) 

def readFile_calculateAF(infile):
    f = open(infile)
    n = len(f.readline().strip().split())
    f.close()
    f = open(infile)
    af = [sum(map(float,line.strip().split()))/(2*n) for line in f.readlines()]
    f.close()
    return(af)

def readFile_process(infile):
    dic = {}
    f = open(infile)
    for line in f.readlines():
        tmp = line.strip().split('\t')
        dic['Chr'+tmp[0]+':'+tmp[1]] = ''
    return(dic)

def calculate_maf(snpid, af0):
    af1 = [1 - i for i in af0]
    maf = [min(af0[i], af1[i]) for i in range(len(af0))]
    gf0 = [i ** 2 for i in af0]
    gf1 = [2 * af0[i] * af1[i] for i in range(len(af0))]
    gf2 = [i ** 2 for i in af1]
    dic = dict(zip(snpid, [list(a) for a in zip(af0, af1, maf, gf0, gf1, gf2)]))
    return(dic)

def combine_dicts(dict_snp_af, dict_snp_en, db):
    for key in dict_snp_af.keys():
        if key in dict_snp_en.keys():
            dict_snp_af.setdefault(key, []).append('TRUE')
        else:
            dict_snp_af.setdefault(key, []).append('FALSE')
    return(dict_snp_af)

def write_file(content, outfile):
    fw = open(outfile, 'w')
    for key, value in content.items():
        fw.write(key + '\t' + '\t'.join(map(str, value)) + '\n')
    fw.close()
    return None

def result():
    af       = readFile_calculateAF(INFILE1)
    snpid    = read_file(INFILE2)
    dict_snp_en = readFile_process(INFILE3)
    db = INFILE3.split('.')[1]
    dict_snp_af = calculate_maf(snpid, af)
    result = combine_dicts(dict_snp_af, dict_snp_en, db)
    write_file(result, OUTFILE)

if __name__ == '__main__':

    result()




