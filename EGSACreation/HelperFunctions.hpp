// Copyright [2018] Isla Carson
#ifndef MAIN_HELPERFUNCTIONS_HPP_
#define MAIN_HELPERFUNCTIONS_HPP_


#include <deque>
#include <iostream>
#include <fstream>
#include <string>

inline void split_reads(std::deque<std::string> &read_frags,
                        const std::string &read, int depth);
inline double get_percent(const double numerator,
                          const uint_least64_t demoninator);
inline bool file_exists(const std::string &str);
inline bool file_can_exist(const std::string &str);
inline int get_lines_per_read(char first_char);


int get_lines_per_read(char first_char) {
    int lines_per_read;

    if (first_char == '@') {
        // File is in fasta format
        lines_per_read = 4;
    } else if (first_char == '>') {
        // File is in fastq format.
        lines_per_read = 2;
    } else {
        std::cerr << "Fastx file cannot be processed."
                  << " Ensure the first char is either '@' or '>'."
                  << std::endl;
        exit(1);
    }
    return lines_per_read;
}


double get_percent(const double numerator, const uint_least64_t demoninator) {
    double result = numerator / static_cast<double>(demoninator);
    result *= 100;
    return result;
}

bool file_exists(const std::string &str) {
    std::ifstream fs(str);
    bool result = fs.is_open();
    fs.close();
    return result;
}

bool file_can_exist(const std::string &str) {
    std::ofstream fs;
    fs.open(str, std::ofstream::out);
    bool result = fs.is_open();
    fs.close();
    return result;
}

void split_reads(std::deque<std::string> &read_frags,
                 const std::string &read, int depth) {
    std::string current_frag = "";
    for (auto &base : read) {
        if (base != 'N') {
            // Append to current fragment
            current_frag += base;
        } else {
            // end current fragment, and queue if long enough.
            if (current_frag.length() >= depth) {
                // Add termination char
                current_frag += '.';
                read_frags.push_back(current_frag);
            }
            current_frag.clear();
        }
    }

    // Test last remaining fragment for length
    if (current_frag.length() >= depth) {
        //  Add termination char
        current_frag += '.';
        read_frags.push_back(current_frag);
    }
}

#endif  // MAIN_HELPERFUNCTIONS_HPP_
