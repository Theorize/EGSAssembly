// Copyright [2018] Isla Carson
#ifndef ASSEMBLER_BWT_HPP_
#define ASSEMBLER_BWT_HPP_

#include <stdint.h>
#include <bitset>
#include <vector>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>

std::string reverse_complement(std::string *forward);
char complement(char base);

struct Suffix_Interval {
    int_least64_t lower;
    int_least64_t upper;
    int match_length;

    int_least64_t size() const;
    bool valid_size(const int &min) const;
    bool valid() const;

    Suffix_Interval() {
        lower = 0;
        upper = 0;
        match_length = 0;
    }

    Suffix_Interval(int_least64_t l, int_least64_t u) {
        lower = l;
        upper = u;
        match_length = 0;
    }

    Suffix_Interval(int_least64_t l, int_least64_t u, int ml) {
        lower = l;
        upper = u;
        match_length = ml;
    }
};

inline int_least64_t Suffix_Interval::size() const {
    return  0 > upper-lower+1 ? 0 :  upper-lower+1;
}

inline bool Suffix_Interval::valid_size(const int &min) const {
    return upper-lower+1 >= min;
}

inline bool Suffix_Interval::valid() const {
    return lower <= upper;
}

inline std::string reverse_complement(const std::string &forward) {
    int length = forward.length();
    std::string rc(length, '.');
    for (int i = 0; i < length; i++) {
        rc[length-1-i] = complement(forward[i]);
    }
    return rc;
}


inline char complement(char base) {
    if (base == 'A') { return 'T';
    } else if (base == 'C') {return 'G';
    } else if (base == 'G') {return 'C';
    } else if (base == 'T') {return 'A';
    } else if (base == 'a') {return 't';
    } else if (base == 'c') {return 'g';
    } else if (base == 'g') {return 'c';
    } else if (base == 't') {return 'a'; }

    return base;
}



class BWT {
 private:
    std::vector<int_least64_t> m_total_smaller;
    std::vector<int_least64_t> m_fm_index;
    std::vector<bool> m_lcp_flags;
    std::vector<bool> m_visited;
    int_least64_t m_size;
    int m_sort_depth;
    int_least64_t m_num_fragments;
    int_least64_t m_num_classes;
    int m_min_coverage;
    int m_max_coverage;
    std::string m_results_folder;

 public:
    // Constructors
    BWT() = default;
    explicit BWT(const char* input_file,
                 int min_coverage = 70, int max_coverage = 135,
                 const char* results_folder = "");
    ~BWT();

    Suffix_Interval new_SI() const;

    void read_input_file(const char* path);
    std::vector<int_least64_t> comma_delimit_to_nums(std::string line);
    int_least64_t meta_info_to_num(std::string line);

    int convert_base_to_class(char base) const;
    char convert_class_to_base(int char_eq) const;
    char convert_class_to_complemented_base(int char_eq) const;

    // Getters
    int_least64_t get_length() const;
    int_least64_t get_total_smaller_than_class(const int i) const;
    int_least64_t get_cumulative_count(int_least64_t i, int char_eq) const;
    int get_sort_depth() const;
    bool get_lcp_flag(int_least64_t i) const;
    bool get_visited(int_least64_t i) const;
    int get_min_coverage();
    int get_max_coverage();

    void set_max_min_coverage(int min, int max);

    // String matching functions.
    Suffix_Interval update_backwards(Suffix_Interval si, int char_eq) const;
    Suffix_Interval update_backwards(Suffix_Interval si, char c) const;
    void update_backwards(Suffix_Interval *si, int char_eq) const;

    bool is_valid(const Suffix_Interval &si) const;
    bool contains_read(const std::string word) const;

    Suffix_Interval drop_last_base(Suffix_Interval si);
    void drop_last_base(Suffix_Interval *si);

    int8_t find_next_base(const Suffix_Interval &si);
    int8_t find_next_base(const Suffix_Interval &si,
                          std::ofstream &base_dist_output);

    Suffix_Interval get_SI_for(const std::string &kmer) const;
    Suffix_Interval get_SI_for_seed(const std::bitset<160> &kmer_seed) const;

    Suffix_Interval get_SI_for_seed(const std::bitset<160> &kmer_seed,
                                std::vector<Suffix_Interval> *history,
                                int match_length) const;

    Suffix_Interval get_reverse_complement_SI_for(
                                    const std::bitset<160> &kmer_seed) const;
    std::string find_contig(const std::string &kmer);
    std::string find_contig(Suffix_Interval forward, Suffix_Interval reverse,
                            const std::string &kmer,
                            std::ofstream &base_dist_output,
                            std::ofstream &contig_end_reasons);

    std::vector<int8_t> match_procedure(Suffix_Interval &si);
    std::vector<int8_t> match_procedure(Suffix_Interval &si,
                                        std::ofstream &base_dist_output,
                                        std::string &contig_end_reason);
    void clear_visited_log();

    void generate_all_contigs(int expected_repeat_num = 0);
    void mark_reverse_as_visited(const std::vector<int8_t> &contig,
                                 const std::string &kmer,
                                 bool kmer_needs_reverse_complemented = 0,
                                 int_least64_t i = -1);


    std::string convert_to_kmer(std::bitset<160> kmer_seed) const;
    char get_base_at(int i, const std::bitset<160> &kmer_seed) const;
    int get_char_eq_at(int i, const std::bitset<160> &kmer_seed) const;
    int add_one_at_index(int n, std::bitset<160> *kmer_seed) const;
    bool max_seed_passed(const std::bitset<160> &kmer_seed);
};


#endif  // ASSEMBLER_BWT_HPP_
