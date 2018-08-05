// Copyright [2018] Isla Carson
#include <stdint.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <ctime>

#include "BWT.hpp"



BWT::BWT(const char* input_file,
         int min_coverage, int max_coverage,
         const char* results_folder) {
    // Initialise class members
    m_min_coverage = min_coverage;
    m_max_coverage = max_coverage;

    m_total_smaller.assign(6, 0);

    m_results_folder = results_folder;
    // Check the output files can be opened.
    std::ofstream contigs_output(m_results_folder + "contigs.txt");
    if (!contigs_output.is_open()) {
        std::cerr << "The contigs output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    contigs_output.close();

    std::ofstream base_dist_output(m_results_folder + "base_dist.txt");
    if (!base_dist_output.is_open()) {
        std::cerr << "The contigs output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    base_dist_output << "A,C,G,T,base\n";
    base_dist_output.close();

    std::ofstream length_output(m_results_folder + "contig_length.txt");
    if (!length_output.is_open()) {
        std::cerr << "The contigs output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    length_output.close();

    std::ofstream contig_end_output(m_results_folder
                                            + "contig_end_reasons.txt");
    if (!contig_end_output.is_open()) {
        std::cerr << "The contigs output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    contig_end_output.close();

    // Time the read_input_filefn
    std::clock_t start;
    start = std::clock();
    // Read the input file
    read_input_file(input_file);
    std::cout << "\n----- File read (in "
              << ( std::clock() - start ) / (double) CLOCKS_PER_SEC
              << " secs) -----\n";

    // Calculate the start indice for each char in the SA.
    m_total_smaller[0] = 0;
    for (int i = 1; i < 6; ++i) {
        m_total_smaller[i] = m_total_smaller[i-1]
                                    + m_fm_index[5*(m_size)+(i-1)];
    }

    m_visited.assign(m_size, 0);
}


Suffix_Interval BWT::new_SI() const {
    return {0, m_size-1, 0};
}

void BWT::read_input_file(const char* path) {
    // Open input file and validate it exists.
    std::fstream enhanced_SA_file(path);
    if (!enhanced_SA_file.is_open()) {
        std::cerr << "Enhanced SA file does not exist."
                  << std::endl;
        exit(1);
    }


    std::string line;
    std::getline(enhanced_SA_file, line);  // Skip header line
    std::getline(enhanced_SA_file, line);  // Skip read origin line (not used)

    std::getline(enhanced_SA_file, line);
    m_sort_depth = meta_info_to_num(line);
    std::getline(enhanced_SA_file, line);
    m_size = meta_info_to_num(line);
    std::getline(enhanced_SA_file, line);
    m_num_fragments = meta_info_to_num(line);
    std::getline(enhanced_SA_file, line);
    m_num_classes = meta_info_to_num(line);

    m_lcp_flags.assign(m_size+1, 0);
    m_fm_index.assign(5*(m_size+1), 0);

    std::getline(enhanced_SA_file, line);  // Skip seperation line
    std::getline(enhanced_SA_file, line);  // Skip header line

    int_least64_t line_counter = 0;
    while (std::getline(enhanced_SA_file, line) && line != "\n") {
        // Split string by commas
        std::vector<int_least64_t> numbers = comma_delimit_to_nums(line);
        m_lcp_flags[line_counter]= numbers[0];
        for (int i = 0; i < 5; i++) {
            m_fm_index[line_counter*5 + i] = numbers[i+1];
        }
        line_counter++;
        if (line_counter%1000000 == 0) {
            std::cout << " "
                      << static_cast<double>(line_counter)
                                /static_cast<double>(m_size) * 100
                      << "\%" << std::flush; }
    }
    enhanced_SA_file.close();
    m_size = line_counter - 1;
}


void BWT::set_max_min_coverage(int min, int max) {
    m_min_coverage = min;
    m_max_coverage = max;
}

std::vector<int_least64_t> BWT::comma_delimit_to_nums(std::string line) {
    std::vector<int_least64_t> comma_delimited_line;

    std::string current_number = "";
    for (auto &num : line) {
        if (num != ',') {
            // Append to current fragment
            current_number += num;
        } else {
            // Add current number to the queue.
            comma_delimited_line.push_back(std::stoull(current_number));
            current_number = "";
        }
    }
    comma_delimited_line.push_back(std::stoull(current_number));
    return comma_delimited_line;
}

int_least64_t BWT::meta_info_to_num(std::string line) {
    std::string current_number = "";
    for (auto &num : line) {
        if (num != '-') {
            // Append to current fragment
            current_number += num;
        } else { break; }
    }
    return std::stoull(current_number);
}

BWT::~BWT() {
}

int_least64_t BWT::get_length() const {
    return m_size;
}

int_least64_t BWT::get_total_smaller_than_class(int i) const {
    return m_total_smaller[i];
}


int_least64_t BWT::get_cumulative_count(int_least64_t i, int char_eq) const {
    int_least64_t test = m_fm_index.at(5*i+static_cast<int_least64_t>(char_eq));
    return test;
}

bool BWT::get_visited(int_least64_t i) const {
    return m_visited[i];
}

int BWT::get_sort_depth() const {
    return m_sort_depth;
}

bool BWT::get_lcp_flag(int_least64_t i) const {
    return m_lcp_flags[i];
}

Suffix_Interval BWT::update_backwards(Suffix_Interval si, int char_eq) const {
    /* Update the suffix interval to match char_eq + current suffix */
    // Calc the next SI.
    update_backwards(&si, char_eq);
    return si;
}

void BWT::update_backwards(Suffix_Interval *si, int char_eq) const {
    /* Update the suffix interval to match char_eq + current suffix */
    // Calc the next SI.
    si->lower = m_total_smaller[char_eq] + m_fm_index[5*si->lower+char_eq];
    si->upper = m_total_smaller[char_eq] +
                                    m_fm_index[5*(si->upper+1)+char_eq] - 1;
    si->match_length++;
}

Suffix_Interval BWT::update_backwards(Suffix_Interval si, char c) const {
    /* Update the suffix interval to match c + current suffix */
    int char_eq = convert_base_to_class(c);
    update_backwards(&si, char_eq);
    return si;
}

Suffix_Interval BWT::drop_last_base(Suffix_Interval si) {
    drop_last_base(&si);
    return si;
}

void BWT::drop_last_base(Suffix_Interval *si) {
    /* Update the SI to match the current suffix without the last base. */
    bool curr_lcp_flag = m_lcp_flags[si->lower];
    int_least64_t new_lower = si->lower - 1;
    // Scan upwards through the LCP flags to find the new lower limit of the SI
    while (new_lower >=0 && m_lcp_flags[new_lower] == curr_lcp_flag)
        { new_lower--; }

    int_least64_t new_upper = si->upper + 1;
    while (new_upper <= get_length()
            && m_lcp_flags[new_upper] == curr_lcp_flag)
        { new_upper++; }

    si->lower = new_lower + 1;  // While loop goes one too far.
    si->upper = new_upper - 1;  // While loop goes one too gitfar.
    si->match_length--;
}

Suffix_Interval BWT::get_SI_for(const std::string &kmer) const {
     /*  Find the SI of the read */
    Suffix_Interval current_si = new_SI();
    int i = kmer.length() - 1;

    // Backwards search through the kmer, terminating if not found.
    while (current_si.valid() && i >= 0) {
        current_si = update_backwards(current_si, kmer[i]);
        i -= 1;
    }
    return current_si;
}

Suffix_Interval BWT::get_SI_for_seed(const std::bitset<160> &kmer_seed) const {
     /*  Find the SI of the read */
    Suffix_Interval current_si = new_SI();
    int i = m_sort_depth - 1;
    // Backwards search through the kmer (ie froward through the bitset)
    // , terminating if not found.
    // Note: SIGNIFICANTLY slower testing if valid_size(m_min_coverage) used.
    while (current_si.valid()
            && i >= 0) {
        current_si = update_backwards(current_si,
                                      get_char_eq_at(i, kmer_seed));
        i -= 1;
    }
    return current_si;
}


Suffix_Interval BWT::get_SI_for_seed(const std::bitset<160> &kmer_seed,
                                std::vector<Suffix_Interval> *history,
                                int match_length) const {
    // From the current match length, go forward.
    int i = match_length;
    // Backwards search through the kmer (ie froward through the bitset)
    // , terminating if not found.
    while ((*history)[i + 1].valid()
            && (*history)[i + 1].size() >= m_min_coverage
            && i >= 0) {
        (*history)[i] = update_backwards((*history)[i + 1],
                                      get_char_eq_at(i, kmer_seed));
        i -= 1;
    }
    Suffix_Interval current_si = (*history)[i+1];
    return current_si;
}

Suffix_Interval BWT::get_reverse_complement_SI_for(
                                const std::bitset<160> &kmer_seed) const {
     /*  Find the SI of the reverse complement of the read */
    Suffix_Interval current_si = new_SI();
    int i = 0;

    // Forward search through the complemented kmer (ie backwards through the
    // bitset), terminating if not found.
    while (current_si.valid() && i < m_sort_depth) {
        current_si = update_backwards(current_si,
                                        5 - get_char_eq_at(i, kmer_seed));
        i += 1;
    }
    return current_si;
}


int BWT::get_min_coverage() {
    return m_min_coverage;
}

int BWT::get_max_coverage() {
    return m_max_coverage;
}


bool BWT::contains_read(const std::string read) const {
    /*  Determine if read is contained within the dictionary. */
    Suffix_Interval current_si = new_SI();
    int i = read.length() - 1;

    // Backwards search through the read, terminating if not found.
    while (current_si.lower <= current_si.upper && i >= 0) {
        current_si = update_backwards(current_si, read[i]);
        i -= 1;
    }

    // lower <= upper IFF the full read has been found in the dictionary.
    return (current_si.lower <= current_si.upper);
}


bool BWT::is_valid(const Suffix_Interval &si) const {
    // Note: if the SI is not valid, it has size 0 < m_min_coverage.
    return si.size() >= m_min_coverage
           && si.size() <= m_max_coverage
           && !m_visited[si.lower];
}

int BWT::convert_base_to_class(char base) const {
    switch (base) {
            case 'a': return 1;
            case 'c': return 2;
            case 'g': return 3;
            case 't': return 4;
            case '.': return 0;
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'T': return 4;
            default: return -1;
        }
}

char BWT::convert_class_to_base(int char_eq) const {
    switch (char_eq) {
        case 1: return 'a';
        case 2: return 'c';
        case 3: return 'g';
        case 4: return 't';
        case 0: return '.';
        default: return 'x';
    }
}

char BWT::convert_class_to_complemented_base(int char_eq) const {
    switch (char_eq) {
        case 1: return 't';
        case 2: return 'g';
        case 3: return 'c';
        case 4: return 'a';
        case 0: return '.';
        default: return 'x';
    }
}

std::string BWT::find_contig(const std::string &kmer) {
    Suffix_Interval normal_SI = get_SI_for(kmer);
    std::vector<int8_t> first_contig_half = match_procedure(normal_SI);

    std::string rc_kmer = reverse_complement(kmer);
    Suffix_Interval rc_SI = get_SI_for(rc_kmer);
    std::vector<int8_t> second_contig_half = match_procedure(rc_SI);

    // Mark the reverse complement of the discovered contig as visited.
    if (!first_contig_half.empty()) {
        mark_reverse_as_visited(first_contig_half, kmer); }
    if (!second_contig_half.empty()) {
        mark_reverse_as_visited(second_contig_half, kmer, 1); }
    // Create the contig from the results.
    std::string contig(first_contig_half.size() + kmer.length()
                                            + second_contig_half.size(), '.');
    // First half is in the wrong direction.
    for (int_least64_t i = 0; i < first_contig_half.size(); i++) {
        contig[first_contig_half.size()-1-i] =
                                convert_class_to_base(first_contig_half[i]);
    }
    // kmer is in the correct direction.
    for (int_least64_t i = 0; i < kmer.length(); i++) {
        contig[first_contig_half.size() + i] = toupper(kmer[i]);
    }

    // Second half is in the right direction, but is complemented.
    for (int_least64_t i = 0; i < second_contig_half.size(); i++) {
        contig[first_contig_half.size() + kmer.length() + i]
                = convert_class_to_complemented_base(second_contig_half[i]);
    }

    return contig;
}


std::string BWT::find_contig(Suffix_Interval forward,
                             Suffix_Interval  reverse,
                             const std::string &kmer,
                             std::ofstream &base_dist_output,
                             std::ofstream &contig_end_reasons) {
    std::string first_contig_end_cause = "";
    std::string second_contig_end_cause = "";
    std::vector<int8_t> first_contig_half = match_procedure(forward,
                                                    base_dist_output,
                                                    first_contig_end_cause);
    std::vector<int8_t> second_contig_half = match_procedure(reverse,
                                                    base_dist_output,
                                                    second_contig_end_cause);
    if ( first_contig_half.size() + second_contig_half.size() > 0 ) {
        contig_end_reasons << first_contig_end_cause << "\n"
                           << second_contig_end_cause << "\n";
    }

    // Create the contig from the results.
        std::string contig(first_contig_half.size() + kmer.length()
                                            + second_contig_half.size(), '.');

    if (first_contig_half.size() + second_contig_half.size() > 0) {
         // Mark the reverse complement of the discovered contig as visited.
        mark_reverse_as_visited(first_contig_half, kmer, 0);
        mark_reverse_as_visited(second_contig_half, kmer, 1);

        // First half is in the wrong direction.
        for (int_least64_t i = 0; i < first_contig_half.size(); i++) {
            contig[first_contig_half.size()-1-i] =
                                    convert_class_to_base(first_contig_half[i]);
        }
        // kmer is in the correct direction.
        for (int_least64_t i = 0; i < kmer.length(); i++) {
            contig[first_contig_half.size() + i] = toupper(kmer[i]);
        }

        // Second half is in the right direction, but is complemented.
        for (int_least64_t i = 0; i < second_contig_half.size(); i++) {
            contig[first_contig_half.size() + kmer.length() + i]
                    = convert_class_to_complemented_base(second_contig_half[i]);
        }
    } else {
        contig = kmer;
    }

    return contig;
}


int8_t BWT::find_next_base(const Suffix_Interval &si,
                           std::ofstream &base_dist_output) {
    int8_t next_base;
    int_least64_t max_count = 0;

    int_least64_t lower = 5*si.lower;
    int_least64_t upper = 5*(si.upper+1);

    double total_bases = 0;
    for (int8_t char_eq = 1; char_eq < 5; char_eq++) {
        //  if (upper bound - lower bound) { update max-count and next_base }
        int_least64_t base_count = m_fm_index[upper+char_eq]
                                                - m_fm_index[lower+char_eq];
        total_bases += base_count;
        if (base_count > max_count) {
            max_count = base_count;
            next_base = char_eq;
        }
    }

    std::vector<double> percentages(4, 0.0);
    if ( total_bases == 0.0 ) { return 0; }
    bool not_max_but_big = false;
    for (int8_t char_eq = 1; char_eq < 5; char_eq++) {
        int_least64_t num_bases = m_fm_index[upper+char_eq]
                                                - m_fm_index[lower+char_eq];
        percentages[char_eq-1] = static_cast<double>(num_bases)
                                                            /total_bases;
        if ( (char_eq != next_base && percentages[char_eq-1] >= 0.1) ||
             (char_eq == next_base && percentages[char_eq-1] < 0.8) ) {
            next_base = 0;
        }
        base_dist_output << percentages[char_eq-1] << ",";
    }

    char next_char = convert_class_to_base(next_base);

    // Why do this?  P*genome_len = likely that some kmers have an unusually
    // high number of term_chars next.  This stops the section of the SI (of
    // size total_bases) from falling too far below the minimum allowed
    // coverage.
    if ( total_bases < ((9*m_min_coverage)/10) ) {
        next_base = 6; next_char = 'n'; }

    base_dist_output << next_char << "\n";

    return next_base;
}

int8_t BWT::find_next_base(const Suffix_Interval &si) {
    int8_t next_base;
    int_least64_t max_count = 0;
    int_least64_t max_two = 0;

    int_least64_t lower = 5*si.lower;
    int_least64_t upper = 5*(si.upper+1);
    for (int8_t char_eq = 1; char_eq < 5; char_eq++) {
        //  if (upper bound - lower bound) { update max-count and next_base }
        int_least64_t base_count = m_fm_index[upper+char_eq]
                                                - m_fm_index[lower+char_eq];
        if (base_count > max_count) {
            max_count = base_count;
            next_base = char_eq;
        } else if (base_count == max_count) {
            max_two = base_count;
        }
    }
    if (max_count == max_two) { return 0; }
    return next_base;
}


std::vector<int8_t> BWT::match_procedure(Suffix_Interval &si) {
    std::vector<int8_t> contig;

    int counter = 0;
    while (is_valid(si)) {
        m_visited[si.lower] = 1;
        // Decide next base
        int8_t next_base = find_next_base(si);
        // Mark analysed SI as visited.
        if (next_base != 0) {
            // Update backwards to the new base.
            update_backwards(&si, next_base);
            // Add next base to the contig
            contig.push_back(next_base);
            // Remove the last base from the si
            drop_last_base(&si);
            counter++;
        } else {
            break;
        }
    }

    if (si.size() >= m_min_coverage && si.size() <= m_max_coverage)
    { m_visited[si.lower] = 1; }
    return contig;
}

std::vector<int8_t> BWT::match_procedure(Suffix_Interval &si,
                                         std::ofstream &base_dist_output,
                                         std::string &contig_end_cause) {
    std::vector<int8_t> contig;

    bool break_for_base_ambuiquity = false;
    int8_t next_base = 7;
    while (is_valid(si)) {
        m_visited[si.lower] = 1;
        // Decide next base
        next_base = find_next_base(si, base_dist_output);
        if (next_base != 0 && next_base != 6) {
            // Update backwards to the new base.
            update_backwards(&si, next_base);
            // Remove the last base from the si
            drop_last_base(&si);
            // If the new SI is still valid, push to contig.
            if (is_valid(si)) { contig.push_back(next_base); }
        } else {
            // Note: only reached if this is the first time at this SI
            // Determine and record reason for contig end.
            if (next_base == 0) {
                contig_end_cause.append("ambigiuous result, ");
            } else if (next_base == 6) {
                contig_end_cause.append("too many read ends, ");
            }
            break;
        }
    }

    // Determine reason for contig end, and record.
    if (m_visited[si.lower] == 1 && !(next_base == 6 || next_base == 0)) {
        contig_end_cause.append("already visited, "); }
    if (si.size() < m_min_coverage) {
        contig_end_cause.append("too small, ");
    } else if (si.size() > m_max_coverage) {
        contig_end_cause.append("too large, "); }
    contig_end_cause.append("x\n");

    if (si.size() >= m_min_coverage && si.size() <= m_max_coverage)
        { m_visited[si.lower] = 1; }
    return contig;
}


void BWT::clear_visited_log() {
    m_visited.assign(m_size, 0);
}


void BWT::mark_reverse_as_visited(const std::vector<int8_t> &contig,
                                  const std::string &kmer,
                                  bool kmer_needs_reverse_complemented,
                                  int_least64_t i) {
    // Contig is pulled directly from match_procedure, and so is already
    // reversed.
    /*  Find the SI of the first m_sort_depth bases of the rc-ed contig */
    Suffix_Interval current_si = new_SI();

    // if i and j not provided, default them.
    if (i == -1) { i = contig.size() - 1; }
    // Backwards search through the contig, terminating if not found.
    while ( current_si.match_length < m_sort_depth && i >= 0 ) {
        if ( !current_si.valid() ) {
            if (i + current_si.match_length - 1 >= 0) {
                // Start process again with the next possible k-mer
                mark_reverse_as_visited(contig, kmer,
                                        kmer_needs_reverse_complemented,
                                        i + current_si.match_length - 1);
            }
            // Don't continue down this invalid route.
            return;
        }
        current_si = update_backwards(current_si, 5 - contig[i]);
        i--;
    }

    // Otherwise, if the end of the contig has been reached without the window
    // reaching a match length of m_sort_depth, finish the initial matching.

    // If the kmer needs to be RC-ed before use, this corresponds to using the
    // given kmer, RC-ing it and then RC-ing it again, giving the original.  In
    // This case, the complement doesn't need to be taken and iteration should
    //  be backwards, as normal.
    int j = kmer.length() - 1;
    if (kmer_needs_reverse_complemented) {
        // Known: match_length will reach m_sort_depth before j is too small.
        while (current_si.match_length < m_sort_depth) {
            if ( !current_si.valid() ) {
                if (i + current_si.match_length - (m_sort_depth - 1 - j) - 1 >= 0) {
                    // Start process again with the next possible k-mer
                    mark_reverse_as_visited(contig, kmer,
                    kmer_needs_reverse_complemented,
                    i + current_si.match_length - (m_sort_depth - 1 - j) - 1);
                }
                // Don't continue down this invalid route.
                return;
            }
            int next_base = convert_base_to_class(kmer[j]);  // No complement
            current_si = update_backwards(current_si, next_base);
            j--;
        }
    } else {
        // Known: match_length will reach m_sort_depth before j is too big.
        j = 0;
        while (current_si.match_length < m_sort_depth) {
            if ( !current_si.valid() ) {
                // Start process again with the next possible k-mer
                if (i + current_si.match_length - j - 1>= 0) {
                    mark_reverse_as_visited(contig, kmer,
                                        kmer_needs_reverse_complemented,
                                        i + current_si.match_length - j - 1);
                }
                // Don't continue down this invalid route.
                return;
            }
            int next_base = 5 - convert_base_to_class(kmer[j]);
            current_si = update_backwards(current_si, next_base);
            j++;
        }
    }
    // Now we have a window of length m_sort_depth. Mark as visited.
    if (current_si.size() >= m_min_coverage
                        && current_si.size() <= m_max_coverage)
        { m_visited[current_si.lower] = 1; }


    /* Iterate through the remaining bases in the contig & kmer. */

    // If there's items remaining in the contig, continue matching whilst
    // dropping the end of the window & marking everything as visited.
    while (i >= 0) {
        if ( !current_si.valid() ) {
            if (i + current_si.match_length - 1>= 0) {
            // Start process again with the next possible k-mer
                mark_reverse_as_visited(contig, kmer,
                                        kmer_needs_reverse_complemented,
                                        i + current_si.match_length - 1);
            }
            // Don't continue down this invalid route.
            return;
        }
        // Update backwards to the next base.
        update_backwards(&current_si, 5 - contig[i]);
        // Remove the last base from the si
        drop_last_base(&current_si);
        // Mark analysed SI as visited.
        if (current_si.size() >= m_min_coverage
                        && current_si.size() <= m_max_coverage)
            { m_visited[current_si.lower] = 1; }
        // Reduce i.
        i--;
    }

    // If there's items remaining in the kmer, continue matching whilst
    // dropping the end of the window & marking everything as visited.
    // Note: Our current position in the kmer has already been set above, and
    // the RC of the kmer itself has already been visited.  So j > 0 needed.
    if (kmer_needs_reverse_complemented) {
        while (j > 0) {
            if ( !current_si.valid() ) {
                if (i + current_si.match_length - (m_sort_depth - 1 - j) - 1) {
                // Start process again with the next possible k-mer
                mark_reverse_as_visited(contig, kmer,
                    kmer_needs_reverse_complemented,
                    i + current_si.match_length - (m_sort_depth - 1 - j) - 1);
                }
                // Don't continue down this invalid route.
                return;
            }

            int next_base = convert_base_to_class(kmer[j]);  // No complement
            // Update backwards to the next base.
            update_backwards(&current_si, next_base);
            // Remove the last base from the si
            drop_last_base(&current_si);
            // Mark analysed SI as visited.
            if (current_si.size() >= m_min_coverage
                                && current_si.size() <= m_max_coverage)
                { m_visited[current_si.lower] = 1; }
            // Reduce j.
            j--;
        }
    } else {
        while (j < kmer.length()-1) {
            if ( !current_si.valid() ) {
                if (i + current_si.match_length - j - 1) {
                // Start process again with the next possible k-mer
                mark_reverse_as_visited(contig, kmer,
                                        kmer_needs_reverse_complemented,
                                        i + current_si.match_length - j - 1);
                }
                // Don't continue down this invalid route.
                return;

            }
            int next_base = 5 - convert_base_to_class(kmer[j]);
            // Update backwards to the next base.
            update_backwards(&current_si, next_base);
            // Remove the last base from the si
            drop_last_base(&current_si);
            // Mark analysed SI as visited.
            if (current_si.size() >= m_min_coverage
                                && current_si.size() <= m_max_coverage)
                { m_visited[current_si.lower] = 1; }
            // Reduce j.
            j++;
        }
    }
}



std::string BWT::convert_to_kmer(std::bitset<160> kmer_seed) const {
    std::string result(m_sort_depth, '.');
    for (int i = 0; i < m_sort_depth; i++) {
        result[i] = get_base_at(i, kmer_seed);
    }
    return result;
}

int BWT::add_one_at_index(int n, std::bitset<160> *kmer_seed) const {
    n = 2*n;  // 2 bits per base.

    // Clear all bits below n, to reset to 'a'.
    for (int i = 0; i < n; i++) {
        kmer_seed->reset(i);
    }

    // Flip until first 0.
    int i = n;
    while ( (*kmer_seed)[i] && i < kmer_seed->size() ) {
        kmer_seed->flip(i);
        i++;
    }

    // Flip the rightmost 0 bit.
    kmer_seed->flip(i);

    return i/2;
}

char BWT::get_base_at(int i, const std::bitset<160> &kmer_seed) const {
    int small = 2*i;
    int large = small + 1;
    int8_t ch_eq = 2*kmer_seed[large] + kmer_seed[small] + 1;

    switch (ch_eq) {
        case 1: return 'a';
        case 2: return 'c';
        case 3: return 'g';
        default: return 't';
    }
}

int BWT::get_char_eq_at(int i, const std::bitset<160> &kmer_seed) const {
    int small = 2*i;
    int large = small + 1;
    int result = 2*kmer_seed[large] + kmer_seed[small] + 1;
    return result;
}

void BWT::generate_all_contigs(int expected_repeat_num) {
    // Time the read_input_filefn
    std::clock_t start;
    start = std::clock();

    // Initialise seed:
    std::bitset<160> kmer_seed(0);

    // Open the output files.
    std::ofstream contigs_output(m_results_folder + "contigs.txt",
                                std::ios_base::app);
    if (!contigs_output.is_open()) {
        std::cerr << "The contigs output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    std::ofstream base_dist_output(m_results_folder + "base_dist.txt",
                                std::ios_base::app);
    if (!base_dist_output.is_open()) {
        std::cerr << "The base_dist_output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    std::ofstream length_output(m_results_folder + "contig_length.txt",
                                std::ios_base::app);
    if (!length_output.is_open()) {
        std::cerr << "The length_output file cannot be opened."
                  << std::endl;
        exit(1);
    }
    std::ofstream contig_end_output(m_results_folder
                                            + "contig_end_reasons.txt",
                                    std::ios_base::app);
    if (!contig_end_output.is_open()) {
        std::cerr << "The contig_end_reasons output file cannot be opened."
                  << std::endl;
        exit(1);
    }

    // Set up to provide progress updates.
    uint_least64_t counter = 0;

     // Suffix interval array.
    std::vector<Suffix_Interval> history(m_sort_depth+1);
    for (int i = 0; i < m_sort_depth + 1; i++) {
        history[i] = new_SI();
    }
    int match_length = m_sort_depth - 1;

    Suffix_Interval previous = new_SI();
    while (max_seed_passed(kmer_seed)) {
        Suffix_Interval forward = get_SI_for_seed(kmer_seed, &history,
                                                  match_length);
        if (is_valid(forward)) {
            Suffix_Interval backward = get_reverse_complement_SI_for(kmer_seed);
            std::string kmer = convert_to_kmer(kmer_seed);
            std::string contig = find_contig(forward, backward, kmer,
                                    base_dist_output, contig_end_output);
            if (contig.length() > m_sort_depth) {
                if (expected_repeat_num != 0) {
                    length_output << contig.length() << ","
                                  << expected_repeat_num << "\n";
                } else {
                    length_output << contig.length() << "\n";
                }
                contigs_output << contig << "\n";
            }
        }
        match_length = add_one_at_index(m_sort_depth - forward.match_length,
                                        &kmer_seed);

        counter++;
        if (counter%1000000 == 0) { std::cout << counter/1000000 << " " << std::flush; }
    }

    // Close output files
    contigs_output.close();
    base_dist_output.close();
    length_output.close();
    contig_end_output.close();

    std::cout << "\n----- Assembly completed (in "
          << ( std::clock() - start ) / (double) CLOCKS_PER_SEC
          << " secs, with " << counter << " of 2^"
          << m_sort_depth << " seeds checked) -----\n";
}


bool BWT::max_seed_passed(const std::bitset<160> &kmer_seed) {
    bool can_continue = false || kmer_seed[2*m_sort_depth];
    can_continue = can_continue || kmer_seed[2*m_sort_depth+1];
    can_continue = can_continue || kmer_seed[2*m_sort_depth+2];
    can_continue = can_continue || kmer_seed[2*m_sort_depth+3];
    return !can_continue;
}

