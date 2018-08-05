// Copyright [2018] Isla Carson
#include "LinearSA.hpp"
#include <vector>
#include <bitset>
#include <iterator>
#include <algorithm>
#include <string>
#include <stack>


LinearSA::LinearSA(const std::string &fastx,
                         std::unordered_map<std::string, uint_least64_t> input)
            :ProblemDescription(fastx, input, true) {
    // Process problem and get the equivalence Classes
    m_classes = process_problem_to_vector();
    if ( get_depth() == 0 ) { set_depth_to_max(); }
    is_valid();
    // Character Order
    m_char_order['.'] = 0;
    m_char_order['A'] = 1;
    m_char_order['C'] = 2;
    m_char_order['G'] = 3;
    m_char_order['T'] = 4;

    m_curr_length = 0;
    num_classes = 0;

    // Suffix Array
    m_SA = std::vector<uint_least64_t> (get_SA_length());
}


// Deconstructor
LinearSA::~LinearSA() {}


uint_least64_t LinearSA::get_sort_length() const {
    return m_curr_length;
}

uint_least64_t LinearSA::get_class_at(uint_least64_t i) const {
    return m_classes[i];
}


uint_least64_t LinearSA::get_suffix_at(uint_least64_t i) const {
    return m_SA[i];
}

void LinearSA::initialise_alg() {
    // Iterate through the file by char
    uint_least64_t num_bases = m_char_order.size();
    std::vector<uint_least64_t> count(num_bases, 0);

    for (uint_least64_t i = 0; i < get_SA_length(); i++)
        { count[m_classes[i]] += 1; }

    // Get the starting position for each char in the SA for cyclic shifts
    // of length 1.
    for (uint_least64_t i = 1; i < num_bases; i++) { count[i] += count[i-1]; }

    // Loop through the reads, sorting the first char of each suffix using a
    // count sort
    for (uint_least64_t i = get_SA_length(); i-- > 0;) {
        // Get the class for i.
        uint_least64_t char_standing = m_classes[i];
        // Decrement
        count[char_standing] -= 1;

        // Get it's the suffix location in the SA and assign.
        m_SA[count[char_standing]] = i;
    }

    m_curr_length = 1;
    num_classes = num_bases;
}

void LinearSA::sort_suffix_pairs_count_sort() {
    // Create & initialise the count array.
    // Plus 2 required as counting starts at 0 and we count at the +1 pos.
    std::vector<uint_least64_t> count(get_SA_length(), 0);
    // Create a new vector of the required length.
    std::vector<uint_least64_t> new_SA(get_SA_length(), 0);

    // loop through the full SA, counting the occurence of each letter.
    for (uint_least64_t index = 0; index < get_SA_length(); index += 1) {
        count[m_classes[index]] += 1;
    }

    // Cumulatively total over the key value.
    for (uint_least64_t i = 1; i < get_SA_length(); i++) {
        count[i] += count[i-1];
    }

    for (uint_least64_t i = get_SA_length(); i-- > 0;) {
        uint_least64_t suffix_start = m_SA[i];
        suffix_start = (suffix_start < m_curr_length) ?
                        suffix_start - m_curr_length + get_SA_length() :
                        suffix_start - m_curr_length;

        uint_least64_t suffix_class = m_classes[suffix_start];
        count[suffix_class] -= 1;
        new_SA[count[suffix_class]] = suffix_start;
    }

    // Copy the completed new_SA into m_SA.
    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        m_SA[i] = new_SA[i];
    }
}

void LinearSA::update_classes() {
    std::vector<uint_least64_t> new_classes(get_SA_length());

    // Initialise.
    new_classes[m_SA[0]] = 0;
    uint_least64_t prev_small = m_classes[m_SA[0]];


    uint_least64_t prev_big = m_SA[0] + m_curr_length;
    prev_big = (prev_big >= get_SA_length()) ?
                       m_classes[prev_big - get_SA_length()] :
                       m_classes[prev_big];

    uint_least64_t curr_class = 0;
    for (uint_least64_t i = 1; i < get_SA_length(); i++) {
        uint_least64_t curr_small = m_classes[m_SA[i]];

        uint_least64_t curr_big = m_SA[i] + m_curr_length;
        curr_big = (curr_big >= get_SA_length()) ?
                        m_classes[curr_big - get_SA_length()] :
                        m_classes[curr_big];


        if ( (prev_small != curr_small) || (prev_big != curr_big) ) {
            curr_class += 1;
        }
        new_classes[m_SA[i]] = curr_class;
        prev_small = curr_small;
        prev_big = curr_big;
    }

    for (uint_least64_t j=0; j < get_SA_length(); j++) {
        m_classes[j] = new_classes[j];
    }

    double avg_classes = num_classes;
    num_classes = m_classes[m_SA[get_SA_length()-1]];
    std::cout << "Num new classes: " << num_classes
              << "\nAvg number of unique app_class per prev class:"
              << static_cast<double>(num_classes)/avg_classes
              << "\nAvg number of suffices per current class: "
              << static_cast<double>(get_SA_length())
                                            /static_cast<double>(num_classes)
              << "\n";

    // Increase step length
    m_curr_length *=2;
}


void LinearSA::run_full_alg() {
    // Initialise alg, and first sort from the path.
    initialise_alg();
    std::cout << "\nSearch Completed to depth: \n" <<m_curr_length<< std::endl;

    while (m_curr_length < get_depth()) {
        sort_suffix_pairs_count_sort();
        update_classes();
        std::cout << "Completed to " <<m_curr_length << "\n--------\n"
                  << std::flush;
    }
}

