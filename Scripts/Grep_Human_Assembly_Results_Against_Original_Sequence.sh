st #!/bin/bash
function grep_contigs_30 {
    cat $1/results_for_$2/contigs.txt | while read line
    do
        if [ ${#line} -eq $3 ]
        then
           echo -1 >> $1/results_for_$2/matches.txt
        else
            grep -i -o -m10 $line $1/human30.fa.*source | wc -l >> $1/results_for_$2/matches.txt
        fi
    done
}

function grep_contigs_60 {
    cat $1/results_for_$2/contigs.txt | while read line
    do
        if [ ${#line} -eq $3 ]
        then
           echo -1 >> $1/results_for_$2/matches.txt
        else
            grep -i -o -m10 $line $1/human60.fa.*source | wc -l >> $1/results_for_$2/matches.txt
        fi
    done
}

grep_contigs_60 /isdata/kroghgrp/qpm525/scratch/bac_files/human60 15 16 &
grep_contigs_60 /isdata/kroghgrp/qpm525/scratch/bac_files/human60 35 36 &
grep_contigs_60 /isdata/kroghgrp/qpm525/scratch/bac_files/human60 55 56 &
grep_contigs_60 /isdata/kroghgrp/qpm525/scratch/bac_files/human60 75 76 &
grep_contigs_30 /isdata/kroghgrp/qpm525/scratch/bac_files/human30 15 16 &
grep_contigs_30 /isdata/kroghgrp/qpm525/scratch/bac_files/human30 35 36 &
grep_contigs_30 /isdata/kroghgrp/qpm525/scratch/bac_files/human30 55 56 &
grep_contigs_30 /isdata/kroghgrp/qpm525/scratch/bac_files/human30 75 76 &
wait
 # $1 = genome
 # $2 = kmersize
 # $3 = minimum considered size - k+1 is guaranteed to be in the original seq, so is not checked.
