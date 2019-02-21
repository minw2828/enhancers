#!/usr/bin/env python

################################################################################
#
# Description:
# This module gets the pvalue from GWAS, with the requirement of inputs:
# 1. reference genome annotation file gff.gz
# 2. NCBI reference genome data table which connects two different versions of chromosome name: AC_000158.1 : 01
# 3. geneName list text file, assuming gene name is identical to those in gff.gz file 
# 4. threshold to define significance
#
# Usage:
# python <MODULE_NAME> <infile1> <infile2> <infile3> <threshold> <outfile> 
#
# Author:
# Min Wang (min.wang@depi.vic.gov.au)
#
# Date Created:
# 13 April 2016
# 
# The program is designed to be used with Python 2.6 and 2.7.
#
# Date modified and reason:
#
################################################################################


import sys
import resource
import os
import gzip
import collections


infile1    = sys.argv[1] # genomic.gff.gz
infile2    = sys.argv[2] # NCBI.refgen.infoTable.umd311.tab
infile3    = sys.argv[3] # lesESNPs file 
outfile    = sys.argv[4] # result 


def read_file(infile):
    f = open(infile)
    fc = [line.strip() for line in f.readlines()]
    f.close()
    return(fc)  # no format processing 

def read_gzFile(infile):
    f = gzip.open(infile)
    fc = [line.strip() for line in f.readlines()]
    f.close()
    return(fc) # no format processing

def get_dict_ncbi(fc, sep):
    return(dict([(line.split(sep)[3], line.split(sep)[2]) for line in fc[1:]]))

def filterGffByGeneName(geneName, dict_ncbi, gff):
    tmp = [line.split('\t') for line in gff if 'BestRefSeq' in line and 'gene='+geneName in line] 
    res = []
    for item in tmp:
        if 'AC' in item[0]:
            chr = dict_ncbi[item[0]]
            start = item[3]
            end = item[4]
            gene = [line for line in item[-1].split(';') if 'gene' in line][0].replace('gene=', '')
            description = [i for i in item[-1].split(';') if 'description' in i][0].replace('description=', '')
            res.append([gene, chr, start, end, description])
    return(res) 

def loopOverSigSnps(gff, dict_ncbi, combine):
    res = []
    for item in combine:
        res += filterGffByGeneName(item[0], dict_ncbi, gff)
    return res

def rawProcessLesESNPs(fc):
    return([line.split('\t') for line in fc])

def writeList2File(content, outfile, sep):
    fw = open(outfile, 'w')
    for item in content:
        fw.write(sep.join(map(str, item))+'\n')
    return None

def result():
    gff       = read_gzFile(infile1)
    dict_ncbi = get_dict_ncbi(read_file(infile2), '\t')
    combine  = rawProcessLesESNPs(read_file(infile3))
    res = loopOverSigSnps(gff, dict_ncbi, combine)
    writeList2File(res, outfile, '\t')

if __name__ == '__main__':

    result()




