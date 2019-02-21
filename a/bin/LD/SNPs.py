#!/usr/bin/env python

################################################################################
#
# Description:
# This module calculates pearson r correlation for SNPs in each row.
#
# Usage:
# python <MODULE_NAME> <INPUT_FILE> <OUTPUT_FILE> 
#
# Author:
# Min Wang (min.wang@depi.vic.gov.au)
#
# Date Created:
# 25 April 2016
# 
# The program is designed to be used with Python 2.6 and 2.7.
#
# Date modified and reason:
#
################################################################################


import sys
import os
import itertools
import scipy
from scipy.stats import pearsonr


INFILE     = sys.argv[1]
OUTFILE    = sys.argv[2]


def read_file(infile):
    fc = {}
    f = open(infile)
    for line in f.readlines():
        tmp = line.strip().split(' ')
        fc[tmp[0]] = map(float, tmp[1:])
    f.close()
    return(fc) 

def calculate_pearsonRcorrelation(fc):
    # intepretation of results: http://stats.stackexchange.com/questions/64676/statistical-meaning-of-pearsonr-output-in-python
    Keys = list(itertools.combinations(fc.keys(), 2))
    Values = [scipy.stats.pearsonr(a, b) for a, b in itertools.combinations(fc.values(), 2)]
    return(dict(zip(Keys, Values)))

def writeDict2File(content, outfile):
    fw = open(outfile, 'w')
    for key, value in content.items():
        k = ','.join(map(str, key))
        v = ','.join(map(str, value))
        fw.write(','.join([k, v])+ '\n')
    fw.close()


def result():
    fc = read_file(INFILE)
    res = calculate_pearsonRcorrelation(fc)   
    writeDict2File(res, OUTFILE)

if __name__ == '__main__':

    result()




