// Copyright [2018] Isla Carson
#ifndef MAIN_CONSTRUCTSA_HPP_
#define MAIN_CONSTRUCTSA_HPP_

#include <stdint.h>
#include <map>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <string>
#include <bitset>
#include <mutex>
#include <thread>
#include <stack>
#include "ProblemDescription.hpp"


class ConstructSA: public ProblemDescription {
 private:
    std::vector<uint_least64_t> m_base_indices;
    std::vector<uint_least64_t> m_min_class_for_base;
    int m_alphabet_size;
    std::vector<uint_least64_t> m_SA;
    std::vector<uint_least64_t> m_classes;
    std::vector<bool> m_flags;
    std::vector<bool> m_term_char_passed;
    uint_least64_t m_curr_length;
    uint_least64_t num_classes;

    std::stack<bool> m_alg_steps;        // Alg steps for constructing an SA

    // Limit below which a comparison sort algorithm will be used to sort
    // the radix buckets.
    uint_least64_t m_limit;

    // For threading
    uint_least64_t m_num_threads;
    std::vector<uint_least64_t> m_thread_chunks;
    std::vector<std::thread> m_threads;

    // Setters / finders
    void find_base_indices();

    void set_flag(uint_least64_t i, bool flag);
    void set_sort_choice(uint_least64_t i, bool choice);
    void set_precise_alg_steps();
    void set_full_alg_steps();


    // Do-ers
    void order_and_total(
                std::unordered_map<uint_least64_t, uint_least64_t> &count_map);

    void update_base_class_start();

    uint_least64_t get_classes_in_interval(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                bool (ConstructSA::*get_flag)(uint_least64_t) const,
                uint_least64_t class_start,
                std::vector<uint_least64_t> &class_changes,
                std::vector<bool> &flag_changes);

    uint_least64_t get_immediately_appending_class(uint_least64_t i) const;

    uint_least64_t get_doubled_appending_class(uint_least64_t i) const;

    uint_least64_t get_appending_class(uint_least64_t i, bool mode) const;

    bool get_single_char_flag(uint_least64_t i) const;

    int get_prepending_class(uint_least64_t i) const;

    int get_next_work_chunk(std::vector<bool> *sort_pointer);

    void update_SA_comb_sort(
                uint_least64_t class_start, uint_least64_t class_end,
                std::unordered_map<uint_least64_t, uint_least64_t> &count_map,
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const);

    uint_least64_t get_class_end(uint_least64_t i);

    uint_least64_t count_and_get_class_end(uint_least64_t i,
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                std::unordered_map<uint_least64_t, uint_least64_t> &count_map);

    void sort_thread(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                std::vector<bool> *sort_pointer);

    void get_class_thread(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                bool (ConstructSA::*get_flag)(uint_least64_t) const,
                std::vector<bool> *work_pointer,
                std::vector<std::vector<uint_least64_t>> *class_section,
                std::vector<std::vector<bool>> *flag_section);

    void revise_classes_thread(
                        std::vector<std::vector<uint_least64_t>> *class_section,
                       std::vector<uint_least64_t> *counter,
                       std::vector<bool> *work_pointer);

    void validate_user_input(
                        std::unordered_map<std::string, uint_least64_t> *input);

 public:
    // Constructors
    ConstructSA();

    ConstructSA(const std::string &fastx,
                std::unordered_map<std::string, uint_least64_t> input);

    // Deconstructor
    ~ConstructSA();

    // Getters
    bool get_flag_at(uint_least64_t i) const;

    bool get_sort_choice_for(uint_least64_t i) const;

    uint_least64_t get_base_index(uint_least64_t i) const;

    uint_least64_t get_suffix_at(uint_least64_t i) const;

    uint_least64_t get_class_at(uint_least64_t i) const;

    uint_least64_t get_sort_length() const;

    uint_least64_t get_base_char_at(uint_least64_t i) const;

    uint_least64_t get_min_class_for_base(uint_least64_t i) const;

    std::string see_alg_steps() const;

    std::stack<bool> copy_alg_steps() const;

    // Do-ers
    void initialise_alg();

    void sort_suffix_pairs(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const);

    void update_classes(
        uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
        bool (ConstructSA::*get_flag)(uint_least64_t) const,
        bool mode);

    void run_full_alg();


    // Helper functions for tests.cpp.
    void sort_suffix_pairs(bool mode);

    void update_classes(bool mode);

    void set_all_sort_choices(bool choice);

    void create_flagged_fm_index(std::vector<bool> *lcp_flags,
                                 std::vector<uint_least64_t> *fm_index);

    void create_fm_index(std::vector<uint_least64_t> *fm_index);

    void print_flagged_fm_index(const std::string SA_file);
    void print_fm_index(const std::string SA_file);

    uint_least64_t get_compressed_size();
    std::string create_BWT();
    int convert_base_to_class(char base) const;
    char convert_class_to_base(int char_eq) const;


};



#endif  //  MAIN_CONSTRUCTSA_HPP_
