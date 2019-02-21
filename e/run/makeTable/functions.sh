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

function remove_dir () {
    if [ -d $1 ]; then rm -rf $1; fi
}

function remove_file () {
    if [ -f $1 ]; then rm -rf $1; fi
}

function decompressIfile_redirectOfile () {
    # $1: input .gz file 
    # $2: output directory for decompressed file
    # $3: whether decompress file or not: 0: no; 1: yes, decompress  
    fn_decom=`echo $1 | awk -F'/' '{print $NF}' | sed 's/.gz$//g'`
    if [ $3 -eq 0 ]; then
        :
    elif [ $3 -eq 1 ]; then
        gunzip -c $1 > $2/$fn_decom
    fi
    echo $fn_decom
}

function copyBack () {
    # $1: $target to be copied
    # $2: directory to be copied to
    for i in $1; do if [ -f $i ]; then cp $i $2/; fi; done
    rm -rf $NODEDIR
}

function sortFile_arbitraryOrder () {
    # $1: file to be sorted
    # $2: file with targeted order in the first row
    head -1 $2 | sed 's/\s/\n/g' | awk 'FNR == NR { lineno[$1] = NR; next} {print lineno[$1], $0;}' - $1 | sort -k 1,1n | cut -d' ' -f2- > tmp.pbs
    mv tmp.pbs $1
}

function bool_identical_ids () {
    # $1: an array of filenames. The first columns in these files are compared for identical. 
    target=$1
    for file1 in ${target[@]}; do 
        array1=`awk '{print $1}' $file1`
        for file2 in ${target[@]}; do 
            array2=`awk '{print $1}' $file2`
            if [[ $array1 != $array2 ]]; then 
                printf "IDs in the following two files are not identical.\n$file1\n$file2\n"
                exit 0
            fi 
        done
    done
}

function transpose_file () {
    awk '
    { 
        for (i=1; i<=NF; i++)  {
            a[NR,i] = $i
        }
    }
    NF>p { p = NF }
    END {    
        for(j=1; j<=p; j++) {
            str=a[1,j]
            for(i=2; i<=NR; i++){
                str=str" "a[i,j];
            }
            print str
        }
    }' $1 > $2
}
function transpose_file1 () {
    python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $1 > $2
}

function examine_outDir () {
    # $1: overall output directory 
    # $2: filename of interest 
    for a in $1/[b-c]*; do
        if [ -d $a ]; then
            for b in $a/*; do
                if [ -d $b ]; then
                    for c in $b/*; do
                        if [ -d $c ]; then
                            ls -l $c/$2
                        fi
                    done
                fi
            done
        fi
    done
}

function delete_piled_qsubJobs () {
    # $1: start node number, inclusive
    # $2: end node number, inclusive
    qdel `seq -f "%.0f" $1 $2`
}


