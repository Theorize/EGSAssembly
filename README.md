## Compilation
`cd` into the Assembly or EGSACreation folders and run:
* `make main`  to create the ./main executable
* `make clean` to removes all .o and executable files
* `make tests`:  to create the ./tests executable.  This exe will not run, as test data is not provided.

`cd` into the Scripts folder and run 
* `g++ -std+c++11 -o sim_reads simulate_reads.cpp` then `./sim_reads --help`
* `g++ -std+c++11 -o reverse_complement reverse_complement_fasta.cpp` then `./reverse_complement --help`

## EGSACreation

Setup: `cd ./EGSACreation`, then `make main clean`

Usage: ./main [--help] [--verbose] [command] [parameter=value]

Required: A command and a `-fastx_file` is required.

#### Allowed Commands:
* `probdesc`: Evaluates the fastx file.
* `linear`: Constructs the SA using the linear algorithm
* `threaded`: Constructs the SA using an algorithm suitable for threading.
              
#### Allowed Parameters:
* `-fastx_file`: file path of fastx file to be processed
* `-sort_depth`: Character depth to which the SA will be sorted
* `-min_frag_len`: Minimum acceptable read fragment length
* `-precise_depth`: 0 if SA sorting depth can exceed `-sort_depth`, 1 otherwise.
* `-num_threads`: Number of threads to be used. Min 1, Max 1024. Int only.
* `-alg_limit`: Bucket size below which the threaded alg will switch to comparison sort.  Min 1.  Int only.


#### Defaults
* `-num_threads`: 1
* `-alg_limit`: 80000
1. If `-sort_depth` given:
  * `-min_frag_len` will default to `-sort_depth`.
  * `precise_depth` will default to 1.
2. If `-sort_depth` is not given:
  * `-min_frag_len` will default to 1.
  * `-sort_depth`  will be set to the length of the longest read fragment, so that full sorting is performed.
  * `precise_depth` will default to 1, as sorting to a character length beyond that of the longest read fragment will have no impact on the result.
 
#### Notes
1. `-precise_depth`, `-alg_limit`, and `-num_threads` are only used by `threaded`.
2. The linear algorithm always sorts to a depth 2^n where n is the minimum value such that 2^n >= `-sort_depth`.


## Assembly

Setup: `cd ./Assembly`, then  `make main clean`, then

Useage: ./main [ESGA file] [minimum copy number] [maximum copy number] [results folder location]

Required: An EGSA file created by EGSACreation, with inputs  `-sort_depth` ,   `-precise_depth=1`,    `-include_rc=1`, and   `-output_file`.    
