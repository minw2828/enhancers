#!/usr/bin/env python

################################################################################
#
# Description:
# This module calculates allele frequency and genotypic frequency.
#
# Usage:
# python <MODULE_NAME> <INFILE1> <INFILE2> <OUTFILE>
#
# Author:
# Min Wang (min.wang@depi.vic.gov.au)
#
# Date Created:
# 06 May 2015
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
OUTFILE    = sys.argv[3]


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

def format_output(contentList, sep):
    return(sep.join(map(str,contentList)))   

def writeDict2File(content, outfile, sep):
    fw = open(outfile, 'w')
    for key, value in content.items():
        fw.write(format_output([key, format_output(value, sep)], sep) + '\n')
    fw.close()
    return None

def result():
    af       = readFile_calculateAF(INFILE1)
    snpid    = read_file(INFILE2)
    dict_snp_af = calculate_maf(snpid, af)
    writeDict2File(dict_snp_af, OUTFILE, '\t')

if __name__ == '__main__':

    result()




