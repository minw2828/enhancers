#!/bin/bash

# Description:
# This script checks if the jobs are properly finished.
# Min Wang, Feb 2017

# Requirement:
# You will need to adjust parameters under [program], [input] and [output] to 
# fit your case.

# Execution:
# sh call.enrichment.scr [script_name] [parallel_index]


#################### parameters ####################
projdate='2016-01-19'
psf='a'
task='TFBS'
jobIdPattern=`echo $1 | sed 's/\.$//g'`

#################### Define paths ####################
path_pre='/group/dairy/Min/geno2pheno'
path_analyse="$path_pre/analyses/$projdate/$psf"
path_data="$path_pre/data"
path_software="$path_pre/software"
binpath="$path_analyse/bin"; runpath="$path_analyse/run"; outpath="$path_analyse/out"

#################### analysis ####################
for i in $jobIdPattern.e*; do head -1 $i; done

echo ""

for i in $jobIdPattern.e*; do tail -1 $i; done

echo ""

for i in $jobIdPattern.o*; do grep "=" "$i"; echo ""; done | awk '{print "#", $0}'

exit 0

