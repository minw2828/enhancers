#!/usr/bin/env python

################################################################################
#
# Description:
# This module transponses an input matrix. 
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
# 06 April 2016
# 
# The program is designed to be used with Python 2.6 and 2.7.
#
# Date modified and reason:
#
################################################################################


import os.path
import sys
import time
import resource
import os


INPUT_FILE     = sys.argv[1]
OUTPUT_FILE    = sys.argv[2]


def read_file(infile):
    f = open(infile)
    fc = [line.strip() for line in f.readlines()]
    f.close()
    return(fc) 

def transpose(fc):
    ids = [line.split(' ')[0] for line in fc]
    body = zip(*[line.split(' ')[1:] for line in fc])
    return(ids + body)

def write_file(content, outfile):
    fw = open(outfile, 'w')
    fw.write(' '.join(content[0])+'\n')
    fw.writelines(' '.join(i) + '\n' for i in content[1:])
    fw.close()

def result():
    fc = read_file(INPUT_FILE)
    transponsed_fc = transpose(fc)
    write_file(transponsed_fc, OUTPUT_FILE)

if __name__ == '__main__':

    result()




