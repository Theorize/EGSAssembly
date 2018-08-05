#!/bin/bash

if [[ ! -e ../data/reverse_complement_$1 ]]; then
    echo "here"
    g++ -std=c++11 -Ofast -o rc ./reverse_complement_fasta.cpp
    ./rc ../data/$1 ../data/reverse_complement_$1
    echo "Reverse complement created."
fi

for folder in human30 human60 ERR430993 ERR431029
do
    mkdir ../DataResults/$folder/kmer_stats

    for kmersize in 15 35 55 75
    do
        ../KMC/exes/kmc -k$kmersize -m10 -sm -fa -cs500 -ci5 -t3 ../DataResults/$folder/bac.fasta ../DataResults/$folder/kmer_stats/kmcoutput ../KMC/exes/temp

        ../KMC/exes/kmc -k$kmersize -m10 -sm -fa -cs500 -ci5 -t3 ../DataResults/$folder/bac_rc.fasta ../DataResults/$folder/kmer_stats/kmcoutput_RC ../KMC/exes/temp

        ../KMC/exes/kmc_tools simple ../DataResults/$folder/kmer_stats/kmcoutput ../DataResults/$folder/kmer_stats/kmcoutput_RC union ../DataResults/$folder/kmer_stats/kmcoutput_both

        ../KMC/exes/kmc_tools transform ../DataResults/$folder/kmer_stats/kmcoutput_both histogram ../DataResults/$folder/kmer_stats/hist$kmersize.txt



        rm ../DataResults/$folder/kmer_stats/kmcoutput.kmc_pre
        rm ../DataResults/$folder/kmer_stats/kmcoutput.kmc_suf

        rm ../DataResults/$folder/kmer_stats/kmcoutput_RC.kmc_pre
        rm ../DataResults/$folder/kmer_stats/kmcoutput_RC.kmc_suf
        rm ../DataResults/$folder/kmer_stats/kmcoutput_both.kmc_pre
        rm ../DataResults/$folder/kmer_stats/kmcoutput_both.kmc_suf


    done
done

# ../KMC/exes/kmc -k$kmersize -m10 -sm -fa -ci1 -cs500 -t3 ../DataResults/$1_rc.fasta ../KMC/histograms/koutput_rc$kmersize ../KMC/exes/temp

# ../KMC/exes/kmc_tools simple ../KMC/histograms/kmcoutput ../KMC/histograms/koutput_rc$kmersize union ../KMC/histograms/both_fwd_rc$kmersize

# ../KMC/exes/kmc_tools transform ../KMC/histograms/both_fwd_rc$kmersize histogram ../KMC/histograms/hist$kmersize.txt
