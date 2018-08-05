#!/bin/bash
function diff_forward_contigs_against_reverse {
    cat $1/results_for_$2/contigs.txt | while read line
    do
        if [ $(grep -i -c -m1 $line $1/results_for_$2/reverse_contigs.txt) -eq 0 ]
        then
           echo "$line in contigs but not reverse_contigs"
        fi
    done

    cat $1/results_for_$2/reverse_contigs.txt | while read line
    do
        if [ $(grep -i -c -m1 $line $1/results_for_$2/contigs.txt) -eq 0 ]
        then
           echo "$line in reverse_contigs but not contigs"
        fi
    done
}

diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human60 15 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human60 35 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human60 55 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human60 75 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human30 15 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human30 35 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human30 55 &
diff_forward_contigs_against_reverse /isdata/kroghgrp/qpm525/scratch/bac_files/human30 75 &
wait
 # $1 = genome
 # $1 = kmersize
