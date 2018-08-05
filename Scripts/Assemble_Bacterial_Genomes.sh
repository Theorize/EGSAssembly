root_folder="/binf-isilon/kroghgrp/qpm525/scratch/bac_files"

/usr/bin/time -v -o $root_folder/ERR431029/results_for_15/time.txt ../Assembler/main $root_folder/ERR431029/results_for_15/ESA.fm 112 246 $root_folder/ERR431029/results_for_15/ &

/usr/bin/time -v -o $root_folder/ERR431029/results_for_35/time.txt ../Assembler/main $root_folder/ERR431029/results_for_35/ESA.fm 66 192 $root_folder/ERR431029/results_for_35/ &
wait

/usr/bin/time -v -o $root_folder/ERR431029/results_for_55/time.txt ../Assembler/main $root_folder/ERR431029/results_for_55/ESA.fm 34 136 $root_folder/ERR431029/results_for_55/ &

/usr/bin/time -v -o $root_folder/ERR431029/results_for_75/time.txt ../Assembler/main $root_folder/ERR431029/results_for_75/ESA.fm 10 82 $root_folder/ERR431029/results_for_75/ &
wait

/usr/bin/time -v -o $root_folder/ERR430993/results_for_15/time.txt ../Assembler/main $root_folder/ERR430993/results_for_15/ESA.fm 174 302 $root_folder/ERR430993/results_for_15/ &

/usr/bin/time -v -o $root_folder/ERR430993/results_for_35/time.txt ../Assembler/main $root_folder/ERR430993/results_for_35/ESA.fm 104 240 $root_folder/ERR430993/results_for_35/ &
wait

/usr/bin/time -v -o $root_folder/ERR430993/results_for_55/time.txt ../Assembler/main $root_folder/ERR430993/results_for_55/ESA.fm 56 172 $root_folder/ERR430993/results_for_55/ &

/usr/bin/time -v -o $root_folder/ERR430993/results_for_75/time.txt ../Assembler/main $root_folder/ERR430993/results_for_75/ESA.fm 18 102 $root_folder/ERR430993/results_for_75/ &
wait
