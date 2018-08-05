// Copyright [2018] Isla Carson
#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <queue>
#include <map>
#include <deque>
#include <vector>
#include "ProblemDescription.hpp"
#include "HelperFunctions.hpp"

char complement(char base);
uint_least64_t complement(uint_least64_t base);
std::string get_reverse_complement(std::string holder);


uint_least64_t complement(uint_least64_t base) {
    if (base == 0) { return 0;}
    return 5 - base;
}

char complement(char base) {
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

void get_reverse_complement(std::string *holder) {
    std::reverse(holder->begin(), holder->end());
    for (int i = 0; i < holder->length(); i++) {
        (*holder)[i] = complement((*holder)[i]);
    }
}

// Complete User Input
std::unordered_map<std::string, uint_least64_t>
            ProblemDescription::complete_user_input(
                    std::unordered_map<std::string, uint_least64_t> input) {
    // If depth given, set cutoff to depth unless otherwise specified.
    // 1 if given, 0 otherwise.
    bool depth_specified = (input.find("-sort_depth") != input.end());
    bool cutoff_specified = (input.find("-min_frag_len") != input.end());
    if (!depth_specified && !cutoff_specified) {
        std::cerr
            << "Either sort depth or cutoff must be specified. "
            << "See './main --help'."
            << std::endl;
        exit(1);
    }
    if (depth_specified) {
        // Check depth is >= 4.
        if (input["-sort_depth"] > 500) {
            std::cerr
                << "Sort Depth should be  value should be >= 0 and <= 500 "
                << "or unspecified.  0 signifies maximum sort depth. "
                << "See './main --help'."
                << std::endl;
            exit(1);
        }
        input.emplace("-min_frag_len", input["-sort_depth"]);
        // Precision sorting defaults to true.
        input.emplace("-precise_depth", 1);
    } else {
        input.emplace("-min_frag_len", 1);
        // Signal that -sort_depth is unspecified by setting it to 0.
        input["-sort_depth"] = 0;
        input.emplace("-precise_depth", 0);  // Precision sorting to be false.
        // Check provided value of -precise_depth is not 1.
        // (I.e. Precision sorting is disabled if -sort_depth not given.)
        if (input["-precise_depth"] == 1) {
            std::cerr
                << "If -sort_depth is unspecified, -precise_depth must be 0. "
                << "See './main --help'."
                << std::endl;
            exit(1);
        }
    }

    // Reverse complement requirement defaults to no.
    input.emplace("-include_rc", 0);

    // Check cut off >= 4.
    if (input["-min_frag_len"] < 1) {
        std::cerr << "Min frag length should be > 0 or unspecified. "
                  << "See './main --help'."
                  << std::endl;
        exit(1);
    }

    return input;
}


// Constructors
ProblemDescription::ProblemDescription(const std::string &fastx,
                        std::unordered_map<std::string, uint_least64_t> input,
                        bool from_ConstructSA_class) {
    // Check if the file is a reasonable choice.
    assert(file_exists(fastx));
    // Then assign.
    m_fastx_file = fastx;
    m_read_file = "unassigned";

    input = complete_user_input(input);


    m_min_frag_len = input["-min_frag_len"];
    m_depth = input["-sort_depth"];
    m_precision = input["-precise_depth"];
    m_reverse_complement = input["-include_rc"];

    m_read_number = 0;
    m_total_bases = 0;
    m_useable_bases = 0;
    m_percentage_useful_reads = 0.0;
    m_max_fraq_len = 0;
    m_percentage_useful_bases = 0.0;
    m_avg_frag_length = 0.0;
    m_fragment_number = 0;

    // Process problem if not initialised from Construct_SA constructor.
    if (!from_ConstructSA_class) {
        if (m_depth == 0) {
            m_depth = process_problem();
        } else {
            process_problem();
        }
        bool valid = is_valid();
        assert(valid);
    }
}

ProblemDescription::ProblemDescription(const std::string &fastx,
                        const std::string &read_file,
                        std::unordered_map<std::string, uint_least64_t> input) {
    // Check if the files are reasonable choices.
    assert(file_exists(fastx));
    assert(fastx != read_file);
    // Then assign.
    m_fastx_file = fastx;
    m_read_file = read_file;

    input = complete_user_input(input);

    m_min_frag_len = input["-min_frag_len"];
    m_depth = input["-sort_depth"];
    m_precision = input["-precise_depth"];
    m_reverse_complement = input["-include_rc"];



    m_read_number = 0;
    m_total_bases = 0;
    m_useable_bases = 0;
    m_percentage_useful_reads = 0.0;
    m_max_fraq_len = 0;
    m_percentage_useful_bases = 0.0;
    m_avg_frag_length = 0.0;
    m_fragment_number = 0;
    if (m_depth == 0) {
        m_depth = process_problem_to_file();
    } else {
        process_problem_to_file();
    }
    bool valid = is_valid();
    assert(valid);
}


// Destructor
ProblemDescription::~ProblemDescription() {}

int ProblemDescription::convert_base_to_class(char base) {
    if (base == 'A') { return 1;
    } else if (base == 'C') {return 2;
    } else if (base == 'G') {return 3;
    } else if (base == 'T') {return 4;
    } else if (base == '.') {return 0;
    } else if (base == 'a') {return 1;
    } else if (base == 'c') {return 2;
    } else if (base == 'g') {return 3;
    } else if (base == 't') {return 4; }

    // Incorrect if this far is reached.  Return -1 to signify this.
    return -1;
}

// Setters
void ProblemDescription::set_depths(uint_least64_t depth,
                                    uint_least64_t cutoff) {
    m_depth = depth;
    if (cutoff == 0) {
        m_min_frag_len = depth;
    } else {
        m_min_frag_len = cutoff;
    }
}


// Getters
uint_least64_t ProblemDescription::get_depth() const { return m_depth; }

uint_least64_t ProblemDescription::get_min_frag_len() const {
    return m_min_frag_len;
}

bool ProblemDescription::get_reverse_complement_req() const {
    return m_reverse_complement;
}

uint_least64_t ProblemDescription::get_precision_req() const {
    return m_precision;
}

const std::string ProblemDescription::get_fastx_file_name() const {
    // copy constructor for char *.
    return m_fastx_file;
}

uint_least64_t ProblemDescription::get_num_useable_bases() const {
    return m_useable_bases;
}

uint_least64_t ProblemDescription::get_read_number() const {
    return m_read_number;
}

uint_least64_t ProblemDescription::get_total_bases() const {
    return m_total_bases;
}

double ProblemDescription::get_percentage_useful_reads() const {
    return m_percentage_useful_reads;
}

double ProblemDescription::get_percentage_useful_bases() const {
    return m_percentage_useful_bases;
}

double ProblemDescription::get_avg_frag_length() const {
    return m_avg_frag_length;
}

uint_least64_t ProblemDescription::get_fragment_number() const {
    return m_fragment_number;
}

uint_least64_t ProblemDescription::get_SA_length() const {
    return m_fragment_number + m_useable_bases;
}

uint_least64_t ProblemDescription::get_imprecise_sort_depth() const {
    uint_least64_t result = 1;
    while (result < m_depth) {
        result*=2;
    }
    return result;
}

std::ostream& operator<<(std::ostream& output,
                         const ProblemDescription& prob) {
    std::string file(prob.m_fastx_file);
    // Format information nicely.
    output << "Fastx File: " << file << "\n\nNumber of Reads: "
           << std::to_string(prob.m_read_number)
           << "\nPercentage Useful Reads: "
           << std::to_string(prob.m_percentage_useful_reads)
           << "%\n\nNumber of Fragments: "
           << std::to_string(prob.m_fragment_number)
           << "\nAvg. Fragment Length: "
           << std::to_string(prob.m_avg_frag_length)
           << "\nMax. Fragment Length: "
           << std::to_string(prob.m_max_fraq_len)
           << "\n\nTotal Bases: "
           << std::to_string(prob.m_total_bases)
           << "\nSuffix Array Length: "
           << std::to_string(prob.get_SA_length())
           << "\n\nUseable Bases: "
           << std::to_string(prob.m_useable_bases)
           << "\nPercentage Useful Bases: "
           << std::to_string(prob.m_percentage_useful_bases)
           << "\n\nPrecise Sort Depth: "
           << std::to_string(prob.m_depth)
           << "\nImprecise Sort Depth: "
           << std::to_string(prob.get_imprecise_sort_depth())
           << "\nPrecise Sort Depth Requested: "
           << std::to_string(prob.m_precision)
           << "\nMinimum Fragment Size: "
           << std::to_string(prob.m_min_frag_len)
           << "\nReverse Complement included: "
           << std::to_string(prob.get_reverse_complement_req());
    return output;
}


// Do-ers
uint_least64_t ProblemDescription::process_problem() {
    // Open file and validate it exists.
    std::ifstream fastx(m_fastx_file, std::ios::ate);
    if (!fastx.is_open()) {
        std::cerr << "Fastx file does not exist." << std::endl;
        exit(1);
    }
    size_t length = fastx.tellg();
    fastx.seekg(0, fastx.beg);

    // Determine if file is fasta, fastq or neither  - which exits the program
    uint_least64_t c = fastx.peek();
    uint_least64_t lines_per_read = get_lines_per_read(c);

    // Used to determine if a line is a read or not.
    uint_least64_t line_counter = 0;

    // Read data as a block
    char *contents = new char[length];
    fastx.read(contents, length);
    fastx.close();

    // Check the file size isn't too large.
    // Must have sort_depth additional space to allow for some suffix addition
    // operations
    if (length/(lines_per_read/2) + get_imprecise_sort_depth()
        >= UINT_LEAST64_MAX) {
        std::cerr << "Expected file length is "
                  << length/(lines_per_read/2)
                  << ", which exceeds the maximum int size, "
                  << UINT_LEAST64_MAX
                  << " minus the imprecise sort_depth ("
                  << get_imprecise_sort_depth() << ").";
        exit(1);
    }

    // Process contents, recording stats.
    m_max_fraq_len = 0;
    uint_least64_t file_index = 0;
    // uint_least64_t class_index = 0;
    while (file_index < length) {
        // Reset counters
        uint_least64_t line_index = 0;

        // Get next line
        char* line = contents+file_index;
        if (line_counter%lines_per_read == 1) {
            // This line is a read
            m_read_number += 1;

            bool useful_read = 0;
            while (file_index+line_index < length) {
                // Within the same read
                uint_least64_t frag_index  = 0;
                char* frag = line + line_index;
                while (convert_base_to_class(frag[frag_index]) != -1) {
                    // Within the same read frag.
                    frag_index++;
                }

                // Analyse fragment
                if (frag_index < m_min_frag_len) {
                    // Remove fragment from record.
                } else {
                    // Record stats.
                    useful_read = 1;
                    m_fragment_number += 1;
                    m_avg_frag_length += frag_index;
                    m_useable_bases += frag_index;
                    m_max_fraq_len = (frag_index < m_max_fraq_len)?m_max_fraq_len:frag_index;
                }

                // Skip bad char and continue processing line.
                line_index += frag_index + 1;
                // Check if there's further fragments on the line.
                if (frag[frag_index] == '\n') { break; }
            }

            // Stats for full line length.
            m_percentage_useful_reads += useful_read;
            m_total_bases += line_index-1;

        } else {
            // Not a read, accelerate to line end.
            while (line[line_index] != '\n') { line_index++; }
            line_index++;  // skip '\n' char
        }

        // Move onto the next line.
        line_counter += 1;
        file_index = file_index + line_index;
    }

    // Include the termination char in the max frag len,
    m_max_fraq_len++;

    // Stats changed to %s / avgs
    m_percentage_useful_reads = get_percent(m_percentage_useful_reads,
                                            m_read_number);
    m_percentage_useful_bases = get_percent(m_useable_bases,
                                            m_total_bases);
    m_avg_frag_length /= m_fragment_number;

    if (get_reverse_complement_req()) {
        m_total_bases *= 2;
        m_useable_bases *= 2;
        m_fragment_number *= 2;
    }

    delete[] contents;

    return m_max_fraq_len;
}

std::vector<uint_least64_t> ProblemDescription::process_problem_to_vector() {
    // Open file and validate it exists.
    std::ifstream fastx(m_fastx_file, std::ios::ate);
    if (!fastx.is_open()) {
        std::cerr << "Fastx file does not exist." << std::endl;
        exit(1);
    }
    size_t length = fastx.tellg();
    fastx.seekg(0, fastx.beg);

    // Determine if file is fasta, fastq or neither  - which exits the program
    uint_least64_t c = fastx.peek();
    uint_least64_t lines_per_read = get_lines_per_read(c);

    // Used to determine if a line is a read or not.
    uint_least64_t line_counter = 0;

    // Check the file size isn't too large.
    // Must have sort_depth additional space to allow for some suffix addition
    // operations
    if (length/(lines_per_read/2) + get_imprecise_sort_depth()
        > UINT_LEAST64_MAX) {
        std::cerr << "Expected file length is "
                  << length/(lines_per_read/2)
                  << ", which exceeds the maximum int size, "
                  << UINT_LEAST64_MAX
                  << " minus the imprecise sort_depth ("
                  << get_imprecise_sort_depth() << ").";
        exit(1);
    }

    // Reserve the max used memory.
    // For fasta files this is the length (lines_per_read = 2)
    // For fastq, quality scores are ignored, so we can take len/2.
    std::vector<uint_least64_t> classes;
    classes.assign(length/(lines_per_read/2), 0);

    // Read data as a block
    char *contents = new char[length];
    fastx.read(contents, length);
    fastx.close();

    // Process contents, recording stats.
    m_max_fraq_len = 0;
    uint_least64_t file_index = 0;
    int_least64_t class_index = 0;
    while (file_index < length) {
        // Reset counters
        uint_least64_t line_index = 0;

        // Get next line
        char* line = contents+file_index;
        if (line_counter%lines_per_read == 1) {
            // This line is a read
            m_read_number += 1;

            bool useful_read = 0;
            while (file_index+line_index < length) {
                // Within the same read
                uint_least64_t frag_index  = 0;
                char* frag = line + line_index;
                uint_least64_t letter = convert_base_to_class(frag[0]);
                while (letter != -1) {
                    // Within the same read frag.
                    classes[class_index] = letter;
                    class_index++;
                    frag_index++;
                    letter = convert_base_to_class(frag[frag_index]);
                }

                // Analyse fragment
                if (frag_index < m_min_frag_len) {
                    // Return class_index to correct position
                    class_index = class_index - frag_index;
                } else {
                    // Add termination char
                    classes[class_index] = convert_base_to_class('.');
                    class_index++;

                    // Record stats.
                    useful_read = 1;
                    m_fragment_number += 1;
                    m_avg_frag_length += frag_index;
                    m_useable_bases += frag_index;
                    m_max_fraq_len = (frag_index < m_max_fraq_len)?m_max_fraq_len:frag_index;
                }

                // Skip bad char and continue processing line.
                line_index += frag_index + 1;
                // Check if there's further fragments on the line.
                if (frag[frag_index] == '\n') { break; }
            }

            // Stats for full line length.
            m_percentage_useful_reads += useful_read;
            m_total_bases += line_index-1;

        } else {
            // Not a read, accelerate to line end.
            while (line[line_index] != '\n') { line_index++; }
            line_index++;  // skip '\n' char
        }

        // Move onto the next line.
        line_counter += 1;
        file_index = file_index + line_index;
    }

    // Include the termination char in the max frag len,
    m_max_fraq_len++;

    // Stats changed to %s / avgs
    m_percentage_useful_reads = get_percent(m_percentage_useful_reads,
                                            m_read_number);
    m_percentage_useful_bases = get_percent(m_useable_bases,
                                            m_total_bases);
    m_avg_frag_length /= m_fragment_number;
    delete[] contents;

    if (get_reverse_complement_req()) {
        // Add the rc of each read to the class index.
        classes.resize(2*class_index);

        for (int i = 0; i < class_index; i++) {
            // std::cout << class_index+i << "  " << class_index - i - 1 << "\n" << std::flush;
            classes[class_index+i] = complement(classes[class_index - i - 1]);
        }

        m_total_bases *= 2;
        m_useable_bases *= 2;
        m_fragment_number *= 2;
    } else {
        classes.resize(class_index);
    }


    return classes;
}

uint_least64_t ProblemDescription::process_problem_to_file() {
    // Open read file
    std::ofstream write_output(m_read_file);
    if (!write_output.is_open()) {
        std::cerr << "Read output file cannot be opened." << std::endl;
        exit(1);
    }


    std::ofstream reverse_output(m_read_file + ".rc_temp_file");
    if (!reverse_output.is_open()) {
        std::cerr << "Read temp file cannot be opened." << std::endl;
        exit(1);
    }

    if (!get_reverse_complement_req()) {
        reverse_output.close();
        std::string rc_temp = m_read_file + ".rc_temp_file";
        remove(rc_temp.c_str());
    }

    // Open fastx file and validate it exists.
    std::ifstream fastx(m_fastx_file, std::ios::ate);
    if (!fastx.is_open()) {
        std::cerr << "Fastx file does not exist." << std::endl;
        exit(1);
    }
    size_t length = fastx.tellg();
    fastx.seekg(0, fastx.beg);

    // Determine if file is fasta, fastq or neither  - which exits the program
    char c = fastx.peek();
    uint_least64_t lines_per_read = get_lines_per_read(c);

    // Check the file size isn't too large.
    // Must have sort_depth additional space to allow for some suffix addition
    // operations
    if (length/(lines_per_read/2) + get_imprecise_sort_depth()
        > UINT_LEAST64_MAX) {
        std::cerr << "Expected file length is "
                  << length/(lines_per_read/2)
                  << ", which exceeds the maximum int size, "
                  << UINT_LEAST64_MAX
                  << " minus the imprecise sort_depth ("
                  << get_imprecise_sort_depth() << ").";
        exit(1);
    }

    // Used to determine if a line is a read or not.
    uint_least64_t line_counter = 0;

    // Read data as a block
    char *contents = new char[length];
    fastx.read(contents, length);
    fastx.close();

    // Process contents, recording stats.
    m_max_fraq_len = 0;
    uint_least64_t file_index = 0;
    // uint_least64_t class_index = 0;
    while (file_index < length) {
        // Reset counters
        uint_least64_t line_index = 0;

        // Get next line
        char* line = contents+file_index;
        if (line_counter%lines_per_read == 1) {
            // This line is a read
            m_read_number += 1;

            bool useful_read = 0;
            while (file_index+line_index < length) {
                // Within the same read
                uint_least64_t frag_index  = 0;
                std::string holder = "";
                char* frag = line + line_index;
                while (convert_base_to_class(frag[frag_index]) != -1) {
                    // Within the same read frag.
                    holder += frag[frag_index];
                    frag_index++;
                }

                // Analyse fragment
                if (frag_index < m_min_frag_len) {
                    // Remove fragment from record.
                } else {
                    // Add to file
                    write_output << holder;
                    write_output << '.';
                    if (get_reverse_complement_req()) {
                        get_reverse_complement(&holder);
                        reverse_output << holder;
                        reverse_output << '.';
                    }
                    // Record stats.
                    useful_read = 1;
                    m_fragment_number += 1;
                    m_avg_frag_length += frag_index;
                    m_useable_bases += frag_index;
                    m_max_fraq_len = (frag_index < m_max_fraq_len)?m_max_fraq_len:frag_index;
                }

                // Skip bad char and continue processing line.
                line_index += frag_index + 1;
                // Check if there's further fragments on the line.
                if (frag[frag_index] == '\n') { break; }
            }

            // Stats for full line length.
            m_percentage_useful_reads += useful_read;
            m_total_bases += line_index-1;

        } else {
            // Not a read, accelerate to line end.
            while (line[line_index] != '\n') { line_index++; }
            line_index++;  // skip '\n' char
        }

        // Move onto the next line.
        line_counter += 1;
        file_index = file_index + line_index;
    }

    // Include the termination char in the max frag len,
    m_max_fraq_len++;

    // Stats changed to %s / avgs
    m_percentage_useful_reads = get_percent(m_percentage_useful_reads,
                                            m_read_number);
    m_percentage_useful_bases = get_percent(m_useable_bases,
                                            m_total_bases);
    m_avg_frag_length /= m_fragment_number;

    delete[] contents;

    if (get_reverse_complement_req()) {
        reverse_output.close();
        // Open fastx file and validate it exists.
        std::ifstream reverse_output(m_read_file + ".rc_temp_file", std::ios::ate);
        if (!reverse_output.is_open()) {
            std::cerr << "RC temp file does not exist." << std::endl;
            exit(1);
        }
        length = reverse_output.tellg();
        reverse_output.seekg(0, reverse_output.beg);
        // Read data as a block
        char *reverse_complement_contents = new char[length];
        reverse_output.read(reverse_complement_contents, length);


        reverse_output.close();
        std::string rc_temp = m_read_file + ".rc_temp_file";
        remove(rc_temp.c_str());
        write_output << reverse_complement_contents;
        delete[] reverse_complement_contents;
    }
    // Clean up
    write_output.close();
    return m_max_fraq_len;
}


bool ProblemDescription::is_valid() const {
    /* Problem valid only if m_depth is less than the number of bases.
    *  Note that termination chars are not included in the total base count and
    *  that total bases is 0 if the file hasn't been processed.
    */
    bool valid;
    if (m_fragment_number == 0) {
        std::cerr << "No useable fragments found.  "
                  << "Try increasing min_frag_len."
                  << std::endl;
        exit(1);
    } else if (m_fragment_number == 0 && m_useable_bases < 2) {
        std::cerr << "Only one useable base found."
                  << std::endl;
        exit(1);
    } else if ( m_depth > m_useable_bases ) {
         std::cerr << "Required sort depth is greater than the number of "
                   << "suffices to be sorted."
                   << std::endl;
    }
    valid = m_depth < m_useable_bases;
    valid=valid && !(m_fragment_number == 0 || m_percentage_useful_reads < 10);

    return valid;
}

void ProblemDescription::set_depth_to_max() {
    if (m_max_fraq_len != 0) {
       m_depth = m_max_fraq_len;
    } else {
        std::cerr << "Max fraq length is 0,"
                  << "suggesting the problem has yet to be processed.\n";
        exit(1);
    }
}

void ProblemDescription::reduce_depth_to_max() {
    if (m_max_fraq_len != 0) {
       if (m_depth > m_max_fraq_len) { m_depth = m_max_fraq_len; }
    } else {
        std::cerr << "Max fraq length is 0,"
                  << "suggesting the problem has yet to be processed.\n";
        exit(1);
    }
}

uint_least64_t ProblemDescription::get_max_frag_len() const {
    return m_max_fraq_len;
}

