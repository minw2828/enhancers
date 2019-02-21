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
psf='c'
jobIdPattern=`echo $1 | sed 's/\.$//g'`
if [[ $2 == "" ]]; then 
    errorMsg="Execution halted"
else 
    errorMsg=$2 
fi

#################### Define paths ####################
path_pre='/group/dairy/Min/geno2pheno'
path_analyse="$path_pre/analyses/$projdate/$psf"
path_data="$path_pre/data"
path_software="$path_pre/software"
binpath="$path_analyse/bin"; runpath="$path_analyse/run"; outpath="$path_analyse/out"

#################### analysis ####################
function getHead1 {
    for i in $jobIdPattern.e*; do 
        head -1 $i
    done | uniq 
    echo ""
}

function getTail1 {
    for i in $jobIdPattern.e*; do 
        echo $i
        tail -1 $i
        echo ""
    done
#    echo ""
}

function getTime {
    for i in $jobIdPattern.o*; do 
        grep "=" "$i" 
        echo ""
    done | awk '{print "#", $0}'
    echo ""
}

function getFailOnes () {
    errorMsg=$1
    for i in $jobIdPattern.e*; do
        if grep -q $errorMsg $i; then
            echo $i
            file=`echo $i | sed 's/\.e/\.o/g'`
#            tail -1 $file 
#            echo ""; echo ""; echo ""; echo ""
        fi
    done 
}

#getHead1
#getTail1
#getTime
getFailOnes $errorMsg

exit 0

