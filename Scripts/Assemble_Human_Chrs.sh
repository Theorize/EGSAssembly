# root_folder="/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/DataResults"
root_folder="/binf-isilon/kroghgrp/qpm525/scratch/bac_files"

# # # human30
# /usr/bin/time -v -o $root_folder/human30/results_for_15/time.txt ../Assembler/main $root_folder/human30/results_for_15/ESA.fm 26 27 $root_folder/human30/results_for_15/ human30-15 &
#
# /usr/bin/time -v -o $root_folder/human30/results_for_35/time.txt ../Assembler/main $root_folder/human30/results_for_35/ESA.fm 21 22 $root_folder/human30/results_for_35/ human30-35 &
#
# /usr/bin/time -v -o $root_folder/human30/results_for_55/time.txt ../Assembler/main $root_folder/human30/results_for_55/ESA.fm 16 17 $root_folder/human30/results_for_55/ human30-55 &
#
# /usr/bin/time -v -o $root_folder/human30/results_for_75/time.txt ../Assembler/main $root_folder/human30/results_for_75/ESA.fm 11 12 $root_folder/human30/results_for_75/ human30-75 &
# wait

# # human60
# /usr/bin/time -v -o $root_folder/human60/results_for_15/time.txt ../Assembler/main $root_folder/human60/results_for_15/ESA.fm 53 53 $root_folder/human60/results_for_15/ human60-15 &
#
# /usr/bin/time -v -o $root_folder/human60/results_for_35/time.txt ../Assembler/main $root_folder/human60/results_for_35/ESA.fm 43 43 $root_folder/human60/results_for_35/ human60-35 &
#
# /usr/bin/time -v -o $root_folder/human60/results_for_55/time.txt ../Assembler/main $root_folder/human60/results_for_55/ESA.fm 33 33 $root_folder/human60/results_for_55/ human60-55 &
#
# /usr/bin/time -v -o $root_folder/human60/results_for_75/time.txt ../Assembler/main $root_folder/human60/results_for_75/ESA.fm 23 23 $root_folder/human60/results_for_75/ human60-75 &
# wait



# human30 REVERSE
/usr/bin/time -v -o $root_folder/human30/results_for_15/time.txt ../Assembler/main $root_folder/human30/results_for_15/ESA.fm 26 27 $root_folder/human30/results_for_15/reverse_ human30-15 reverse &

/usr/bin/time -v -o $root_folder/human30/results_for_35/time.txt ../Assembler/main $root_folder/human30/results_for_35/ESA.fm 21 22 $root_folder/human30/results_for_35/reverse_ human30-35 reverse &

/usr/bin/time -v -o $root_folder/human30/results_for_55/time.txt ../Assembler/main $root_folder/human30/results_for_55/ESA.fm 16 17 $root_folder/human30/results_for_55/reverse_ human30-55 reverse &

/usr/bin/time -v -o $root_folder/human30/results_for_75/time.txt ../Assembler/main $root_folder/human30/results_for_75/ESA.fm 11 12 $root_folder/human30/results_for_75/reverse_ human30-75 reverse &
wait

# human60 REVERSE
/usr/bin/time -v -o $root_folder/human60/results_for_15/time.txt ../Assembler/main $root_folder/human60/results_for_15/ESA.fm 53 53 $root_folder/human60/results_for_15/reverse_ human60-15 reverse &

/usr/bin/time -v -o $root_folder/human60/results_for_35/time.txt ../Assembler/main $root_folder/human60/results_for_35/ESA.fm 43 43 $root_folder/human60/results_for_35/reverse_ human60-35 reverse &

/usr/bin/time -v -o $root_folder/human60/results_for_55/time.txt ../Assembler/main $root_folder/human60/results_for_55/ESA.fm 33 33 $root_folder/human60/results_for_55/reverse_ human60-55 reverse &

/usr/bin/time -v -o $root_folder/human60/results_for_75/time.txt ../Assembler/main $root_folder/human60/results_for_75/ESA.fm 23 23 $root_folder/human60/results_for_75/reverse_ human60-75 reverse &
wait
