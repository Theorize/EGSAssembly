// Copyright [2018] Isla Carson
#ifndef MAIN_LINEARSA_HPP_
#define MAIN_LINEARSA_HPP_

#include <map>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <string>
#include "ProblemDescription.hpp"


class LinearSA: public ProblemDescription {
 private:
    std::map<char, uint_least64_t> m_char_order;
    std::vector<uint_least64_t> m_SA;
    std::vector<uint_least64_t> m_classes;
    uint_least64_t m_curr_length;
    uint_least64_t num_classes;

 public:
    // Constructors
    LinearSA();
    LinearSA(const std::string &fastx,
                std::unordered_map<std::string, uint_least64_t> input);

    // Deconstructor
    ~LinearSA();

    // Getters
    uint_least64_t get_suffix_at(uint_least64_t i) const;
    uint_least64_t get_class_at(uint_least64_t i) const;
    uint_least64_t get_sort_length() const;

    // Do-ers
    void initialise_alg();
    void sort_suffix_pairs_count_sort();
    void update_classes();
    void run_full_alg();
};

#endif  // MAIN_LINEARSA_HPP_
