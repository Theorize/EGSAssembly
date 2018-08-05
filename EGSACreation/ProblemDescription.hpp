// Copyright [2018] Isla Carson
#ifndef MAIN_PROBLEMDESCRIPTION_HPP_
#define MAIN_PROBLEMDESCRIPTION_HPP_

#include <stdint.h>

#include <iostream>
#include <stack>
#include <cassert>
#include <string>
#include <queue>
#include <vector>
#include <map>
#include <unordered_map>

class ProblemDescription {
 private:
    std::string m_fastx_file;            // fasta or fastq file
    std::string m_read_file;  // optional / file to write out read frags

    uint_least64_t m_depth;              // sort depth
    uint_least64_t m_min_frag_len;       // lower limit for read frag length
    bool m_precision;                    // sort to precise depth required.
    bool m_reverse_complement;           // true if the RC is to be included.

    uint_least64_t m_read_number;        // Total number of reads in fastx_file
    double m_percentage_useful_reads;    // % reads with useable fragments

    uint_least64_t m_fragment_number;    // Number of read fragments
    double m_avg_frag_length;            // Avg length of the read fragments
    int m_max_fraq_len;

    uint_least64_t m_total_bases;        // Number of bases in fastx_file
    uint_least64_t m_useable_bases;      // Num of bases in the used read frags
    double m_percentage_useful_bases;    // % of total based that are used

    // Setters
    void set_depths(uint_least64_t depth, uint_least64_t cutoff);
    void set_files(const std::string &fastx_file, const std::string &read_file);


 public:
    // Complete User Input
    std::unordered_map<std::string, uint_least64_t> complete_user_input(
                        std::unordered_map<std::string, uint_least64_t> input);

    // Constructors
    ProblemDescription(const std::string &fastx,
                       std::unordered_map<std::string, uint_least64_t> input,
                       bool from_ConstructSA_class = false);
    ProblemDescription(const std::string &fastx, const std::string &read_file,
                       std::unordered_map<std::string, uint_least64_t> input);

    // Destructor
    ~ProblemDescription();

    // Setters
    void set_depth_to_max();
    void reduce_depth_to_max();


    // Getters
    const std::string get_fastx_file_name() const;

    uint_least64_t get_depth() const;
    uint_least64_t get_min_frag_len() const;
    uint_least64_t get_precision_req() const;
    uint_least64_t get_max_frag_len() const;
    bool get_reverse_complement_req() const;

    uint_least64_t get_read_number() const;
    double get_percentage_useful_reads() const;

    uint_least64_t get_fragment_number() const;
    double get_avg_frag_length() const;

    uint_least64_t get_total_bases() const;
    uint_least64_t get_num_useable_bases() const;
    double get_percentage_useful_bases() const;

    uint_least64_t get_SA_length() const;

    // Do-ers (do things)
    int convert_base_to_class(char base);
    uint_least64_t process_problem();  // Returns max fragment len
    std::vector<uint_least64_t> process_problem_to_vector();
    uint_least64_t process_problem_to_file();  // Returns max fragment len
    bool is_valid() const;
    friend std::ostream& operator<<(std::ostream& output,
                                      const ProblemDescription& prob);
    uint_least64_t get_imprecise_sort_depth() const;
};

#endif  // MAIN_PROBLEMDESCRIPTION_HPP_
