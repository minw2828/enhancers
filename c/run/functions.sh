#!/bin/bash

# Description:
# This script supports the other scripts in this directory.
# Min Wang, April 2016

# Execution:
# source [script_name]


#################### functions ####################
function make_dir () { 
    if [ ! -d $1 ]; then mkdir $1; fi 
}

function remove_file () {
    if [ -f $1 ]; then rm -rf $1; fi
}

function sortFile_arbitraryOrder () {
    # $1: file to be sorted
    # $2: file with targeted order in the first row
    head -1 $2 | sed 's/\s/\n/g' | awk 'FNR == NR { lineno[$1] = NR; next} {print lineno[$1], $0;}' - $1 | sort -k 1,1n | cut -d' ' -f2- > tmp.pbs
    mv tmp.pbs $1
}


