#include <stdint.h>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

int get_lines_per_read(char first_char) {
    int lines_per_read;

    if (first_char == '@') {
        // File is in fasta format
        lines_per_read = 4;
    } else if (first_char == '>') {
        // File is in fastq format.
        lines_per_read = 2;
    } else {
        std::cerr << "File processed as bases only."
                  << " Complement taken for every line."
                  << std::endl;
        lines_per_read = 1;
    }
    return lines_per_read;
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

int main(int argc, char *argv[]) {
    std::string input_file;
    std::string output_file_name;

    if (argc == 3) {
        input_file = argv[1];
        output_file_name = argv[2];
    } else {
        std::cerr << "Use format is "
                << "./reverse_complement <input_file> <output_file(optional)>"
                << std::endl;
        exit(0);
    }

    // Open output file
    std::ofstream output_file(output_file_name);
    if (!output_file.is_open()) {
        std::cerr << "Read output file cannot be opened." << std::endl;
        exit(1);
    }

    // Open input file and validate it exists.
    std::ifstream fastx(input_file);
    if (!fastx.is_open()) {
        std::cerr << "Input fastx file does not exist."
                  << std::endl;
        exit(1);
    }

    // Determine if file is fasta, fastq or neither  - which exits the program
    uint_least64_t c = fastx.peek();
    uint_least64_t lines_per_read = get_lines_per_read(c);

    // Used to determine if a line is a read or not.
    uint_least64_t line_counter = 0;

    std::string line;
    while (std::getline(fastx, line)) {

        if (line_counter%lines_per_read != 1 && lines_per_read != 1) {
            // Add to reverse_complement, as is.
            output_file << line + "\n";
        } else {
            for (int i = line.length()-1; i >= 0; i--) {
                output_file << complement(line[i]);
            }
            if (lines_per_read != 1) { output_file << "\n"; }
        }
        line_counter++;
        if (line_counter%100000 == 0) {
            std::cout << line_counter/100000 << ",  " << std::flush;
        }
    }
    std::cout << std::endl;

    fastx.close();
    output_file.close();

}
