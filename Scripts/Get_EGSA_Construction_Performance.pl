#!/bin/bash
# g++ -std=c++11 -o sim_reads simulate_reads.cpp

time_file="construction_time.txt"
results_file="construction_results.txt"
input_file="/isdata/kroghgrp/qpm525/scratch/$1"

# Record the output of the alg to records.txt too.
echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> $time_file
echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >> $results_file
echo -e "File: $1, Threads: 32, Depth: 64, AlgLimit: 80000 \n" >> $time_file
echo -e "File: $1, Threads: 32, Depth: 64, AlgLimit: 80000 \n" >> $results_file
{ /usr/bin/time -v ../main/main --verbose threaded -fastx_file=$input_file -sort_depth=64 -alg_limit=80000 -num_threads=32 -min_frag_len=16 >> $results_file ; } 2>>  $time_file
