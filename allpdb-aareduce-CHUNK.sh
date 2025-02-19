#!/bin/bash
set -u -e

export SEAMLESS_QUEUE_FILE=intermediate/rna-pdbs/.seamless-queue

trap 'kill -1 $(jobs -p); kill $(jobs -p); kill -9 $(jobs -p)' EXIT

listf=intermediate/rna-pdbs/file.list
rm -f $listf-*
j=intermediate/rna-pdbs/jobfile

split -l 200 $listf $listf-
for chunkf in $listf-*; do
    rm -f $SEAMLESS_QUEUE_FILE
    rm -f $j
    echo chunk $chunkf

    seamless-queue &   # start working immediately. You can use -q to get better timings

    for code in $(cat $chunkf); do
        rm -f intermediate/rna-pdbs/$code-aa.pdb.CHECKSUM    
        cmd="aareduce/aareduce.py intermediate/rna-pdbs/$code.pdb --heavy --rna --nalib"
        echo $cmd >> $j
        echo >> $j
        echo >> $j
    done


    echo 'job file generated'
    # You can use -q to get better timings, or -v/-vv to get better progress
    seamless-multi -q --ncores 2 $j

    echo 'Jobs submitted'
    seamless-queue-finish

    rm -f $j
done