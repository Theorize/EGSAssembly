// Copyright [2018] Isla Carson
#include <stdint.h>

#include <iostream>
#include <ctime>
#include <cstring>
#include <set>
#include <fstream>

#include "HelperFunctions.hpp"
#include "LinearSA.hpp"
#include "ConstructSA.hpp"


int process_options(int argc, char *argv[]);
int process_commands(int count, std::string *command, int argc, char *argv[]);
std::unordered_map<std::string, uint_least64_t> process_parameters(int count,
                                                                int argc,
                                                                char *argv[]);
void given_fastx_file(std::string *fastx, std::string *output,
                      int argc, char *argv[]);
void show_manual();

int main(int argc, char *argv[]) {
    // Suppress cout output.
    std::cout.setstate(std::ios_base::failbit);
    // For command
    std::string command;

    // Process command line args.
    int count = process_options(argc, argv);
    count = process_commands(count, &command, argc, argv);
    std::unordered_map<std::string, uint_least64_t> input =
                                        process_parameters(count, argc, argv);
    std::string fastx;
    std::string output_file;
    given_fastx_file(&fastx, &output_file, argc, argv);

    // Depth defaults to 5. (to be deleted later.)
    bool depth_specified = (input.find("-sort_depth") != input.end());
    if (!depth_specified) { input["-sort_depth"] = 5; }

    if (command == "probdesc") {
        std::cout << "----------------------------------------------"
                  << "\nEvaluating problem description and fastx file."
                  << "\n----------------------------------------------"
                  <<  std::endl;
        // Display description of the problem
        ProblemDescription prob(fastx, input);
        std::cout << prob << std::endl;
    } else if (command == "linear") {
        // Create the SA using the linear alg
        std::cout << "---------------------------"
                  << "\nCreating SA via linear alg."
                  << "\n---------------------------"
                  << std::endl;
        LinearSA lin_SA(fastx, input);
        std::cout << lin_SA << std::endl;
        lin_SA.run_full_alg();
    } else if (command == "threaded") {
        std::cout << "-----------------------------"
                  << "\nCreating SA via threaded alg."
                  << "\n-----------------------------"
                  << std::endl;
        // Create the SA using the threaded alg
        ConstructSA Con_SA(fastx, input);
        std::cout << Con_SA << std::endl;
        Con_SA.run_full_alg();

        // Create the FM Index.
        if (output_file != "") {
            if (Con_SA.get_precision_req()) {
                Con_SA.print_flagged_fm_index(output_file);
            } else {
                Con_SA.print_fm_index(output_file);
            }
        }
        Con_SA.get_compressed_size();
    }

    std::cout.clear();
    exit(0);
}


void given_fastx_file(std::string *fastx, std::string *output,
                      int argc, char *argv[]) {
    bool fastx_specified = false;
    bool output_specified = false;

    // Parse command line named parameters.
    for (int count = 1; count < argc; count++) {
        std::string parameter(argv[count]);
        std::size_t delimiter_location = parameter.find_first_of("=");
        if (delimiter_location != std::string::npos) {
            // Split incoming named parameters by the first '=' char.
            std::string before = parameter.substr(0, delimiter_location);
            std::string after = parameter.substr(delimiter_location + 1);

            // Search for the fasta file
            if (before == "-fastx_file") {
                    (*fastx) = after;
                    fastx_specified = true;
            } else if (before == "-output_file") {
                (*output) = after;
                output_specified = true;
            }
        }
    }

    // If fastx unspecified, complain.
    if (!fastx_specified) {
        std::cerr << "Fastx file not specified. "
                  << "See './main --help'."
                  << std::endl;
        exit(1);
    }
    if (!output_specified) { (*output) = ""; }

    if (fastx_specified && !file_exists((*fastx))) {
        std::cerr << "Fastx file cannot be opened."
                  << std::endl;
        exit(1);
    }

    if (output_specified && !file_can_exist((*output))) {
        std::cerr << "Output file cannot be opened."
                  << std::endl;
        exit(1);
    }
}

int process_options(int argc, char *argv[]) {
    if (argc > 12) {
        std::cerr << "Too many command line arguments. "
                  << "See './main --help'."
                  << std::endl;
        exit(1);
    }

    std::set<std::string> allowed_options {"--verbose", "--help"};
    int count = 1;
    while (count < argc && argv[count][0] == '-') {
        // deal with options
        if (std::strcmp(argv[count], "--help") == 0) { show_manual();
        } else if (std::strcmp(argv[count], "--verbose") == 0) {
            std::cout.clear();
        } else {
            std::cerr << "Unknown option: " << argv[count] << "\n. Usage: "
                      << "./main [--help] [--verbose] "
                      << "[command] [<parameter>=<value>]"
                      << std::endl;
            exit(1);
        }
        count++;
    }
    return count;
}

int process_commands(int count, std::string *command, int argc, char *argv[]) {
    std::set<std::string> allowed_commands {"probdesc", "linear", "threaded"};
    // If no command provided, print manual.
    if (count >= argc) {
        std::cerr << "Command must be specified. " << "\nUsage: "
                  << "./main [--help] [--verbose] [command] "
                  << "[<parameter>=<value>]"
                  << std::endl;
        exit(1);
    }

    auto search = allowed_commands.find(argv[count]);
    if (search == allowed_commands.end()) {
        std::cerr << argv[count] << " is not a recognised command. "
                                 << "See './main --help'."
                                 << std::endl;
        exit(1);
    } else {
        command->assign(argv[count]);
    }
    count++;

    // Help allowed for each command.
    if (count < argc) {
        if (std::strcmp(argv[count], "--help") == 0) {
            show_manual();
            // Don't increase count, as either exited, or this could be a param
        }
    }

    return count;
}

std::unordered_map<std::string, uint_least64_t> process_parameters(int count,
                                                                int argc,
                                                                char *argv[]) {
    // Allowed parameter names.
    std::set<std::string> allowed_paras {"-fastx_file",
                             "-num_threads",  "-alg_limit", "-sort_depth",
                             "-precise_depth", "-min_frag_len",
                             "-output_file", "-include_rc"};

    std::unordered_map<std::string, uint_least64_t> input_dict;

    while (count < argc) {
        std::string parameter(argv[count]);
        std::size_t delimiter_location = parameter.find_first_of("=");
        if (delimiter_location == std::string::npos) {
            std::cerr << parameter << " is not a recognised option. "
                      << "See './main --help'."
                      << std::endl;
            exit(1);
        } else {
            // Split incoming named parameters by the first '=' char.
            std::string before = parameter.substr(0, delimiter_location);
            std::string after = parameter.substr(delimiter_location + 1);
            if (after.size() == 0) {
                std::cerr << parameter << " is not specified. "
                          << "See './main --help'."
                          << std::endl;
                exit(1);
            }

            // Ensure the name is in the list of allowed parameters.
            auto search = allowed_paras.find(before);
            if (search == allowed_paras.end()) {
                std::cerr << parameter << " is not a recognised parameter. "
                          << "See './main --help'."
                          << std::endl;
                exit(1);
            } else {
                if (before != "-fastx_file" && before != "-output_file") {
                    // All remaining inputs should be ints
                    // So, check it's not a double.
                    if (after.find_first_not_of("0123456789")
                                                        != std::string::npos) {
                            if (before != "-precise_depth" || before == "-include_rc") {
                                std::cerr << before
                                          << " requires an integer value. "
                                          << "See './main --help'."
                                          << std::endl;
                            } else {
                                std::cerr << before
                                          << " requires either 0 or 1. "
                                          << "See './main --help'."
                                          << std::endl;
                            }
                        exit(1);
                    }
                    // Now try to convert to uint_least64_t.
                    try {
                        if (before == "-precise_depth"
                              || before == "-include_rc") {
                            // If precise depth, should be boolean.
                            if (!(std::stoi(after) == 0
                                    || std::stoi(after) == 1)) {
                                std::cerr << before
                                          << " requires either 0 or 1. "
                                          << "See './main --help'."
                                          << std::endl;
                                exit(1);
                            }
                        }
                        input_dict[before] = std::stoi(after);
                    } catch (...) {
                        std::cerr << before
                                  << " requires an integer. "
                                  << "See './main --help'."
                                  << std::endl;
                        exit(1);
                    }
                }
            }
        }
        count++;
    }
    return input_dict;
}

void show_manual() {
    std::cout.clear();
    std::cout << "Usage: ./main [--help] [--verbose] [command] "
              << "[<parameter>=<value>]"
              << "\n\nHere are the available commands: \n\n"
              << "probdesc         Evaluates the fastx file."
              << "\nlinear           Constructs the SA using a linear algorithm"
              << "\nthreaded         Constructs the SA using an algorithm "
              << "suitable for threading."
              << "\n\nHere are the available parameter and value pairs: \n\n"
              << "-fastx_file      Path to the fasta or fastx file containing\n"
              << "                 the reads to be sorted.\n"
              << "-output_file     Path to the output file."
              << "\n-include_rc    Should the reverse complement of each\n"
              <<   "               read be inlcuded? Bool - 0 (no) or 1 (yes)."
              << "\n-sort_depth      Min depth to which the SA will be sorted\n"
              << "                 Int - Min 1, Max 500, "
              <<                  "Defaults to max fragment length."
              << "\n-min_frag_len    Minimum acceptable fragment length.\n"
              << "                 Int - Min 1, Defaults to sort_depth."
              << "\n-precise_depth   Specifies if a precise sort depth is\n"
              << "                 needed. Bool - 0 (no) or 1 (yes).\n"
              << "                 Default: 0 if sort_depth not specified, 1 "
              <<                  "otherwise."
              << "\n-alg_limit       Max bucket size for using comparison sort"
              << "\n                 Int - Min 1, Defaults to 80000."
              << "\n-num_threads     Number of threads to use.\n"
              << "                 Int - Min 1, Max 1024, Defaults to 1."
              << "\n\nNotes:\n1.  '-fastx_file' is required. "
              << "\n2.  Either '-sort_depth' or '-min_frag_len' must be "
              << "specified."
              << "\n3. '-precise_depth', '-alg_limit', and '-num_threads' are"
              << " only utilised by './main threaded'."
              << std::endl;
    exit(1);
}
