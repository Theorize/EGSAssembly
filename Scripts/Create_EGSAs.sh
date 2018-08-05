#!/bin/bash

function runSAC {
for kmersize in 15 35 55 75
do
    # Make directory for results
    mkdir $1/results_for_$kmersize

    frag_len="$((${kmersize}+1))"
    # Create the fm file
    /usr/bin/time -v ../SAConstruction/main --verbose threaded -fastx_file=$1/$2 -sort_depth=$kmersize -precise_depth=1 -alg_limit=160000 -num_threads=32  -min_frag_len=$frag_len -include_rc=1 2>> ${1}/${kmersize}.txt ;
        echo "$1 $kmersize done"
done
}

runSAC  /binf-isilon/kroghgrp/qpm525/scratch/bac_files/human30 chr1.fasta
runSAC  /binf-isilon/kroghgrp/qpm525/scratch/bac_files/human60 chr1.fasta
runSAC /binf-isilon/kroghgrp/qpm525/scratch/bac_files/ERR430993 bac.fasta
runSAC /binf-isilon/kroghgrp/qpm525/scratch/bac_files/ERR431029 bac.fasta

# $1 = folder
# $2 = filename
