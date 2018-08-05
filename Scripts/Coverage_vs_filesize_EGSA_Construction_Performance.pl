#!/bin/bash

for coverage in 60 40 20
do
    echo "generated/reads_""$coverage""_coverage_""$1""gb.fa"
    perl Get_EGSA_Construction_Performance.pl "generated/reads_""$coverage""_coverage_""$1""gb.fa"
done


# $1 = file size
