#!/bin/bash
g++ -std=c++11 -o sim_reads simulate_reads.cpp

echo "filename,  readlen, coverage, size(gb), bases" > ../data/generated/results.txt

for size in 0.5 1 5 10 15 20
do
    for depth in 10 20 40 60
    do
        echo "$depth" "$size"
        ./sim_reads 120 "$depth" "$size" >> ../data/generated/results.txt
    done
done
