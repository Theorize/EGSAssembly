// Copyright [2018] Isla Carson
#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <string>
#include <vector>
#include <map>

#include "BWT.hpp"
std::map< std::string, std::vector<std::vector<int>> > coverage_regions {
    {"human30-15", {{26, 27}, {52, 54}, {78, 81}, {104, 108},
                    {130, 135}, {156, 162}, {208, 216} }},
    {"human30-35", {{21, 22}, {42, 44}, {63, 66}, {84, 88},
                    {105, 110}, {126, 132}, {168, 176}}},
    {"human30-55", {{16, 17}, {32, 34}, {48, 51}, {64, 68},
                    {80, 85}, {96, 102}}},
    {"human30-75", {{11, 12}, {22, 24}, {33, 36}, {44, 48}, {55, 60}}},

    {"human60-15", {{53-2, 53+2}, {106-2, 106+2}, {159-2, 159+2},
                    {212-2, 212+2}, {265-2, 265+2}, {318-2, 318+2},
                    {371-2, 371+2}, {424-2, 424+2}}},
    {"human60-35", {{43-2, 43+2}, {86-2, 86+2}, {129-2, 129+2}, {172-2, 172+2},
                    {215-2, 215+2}, {258-2, 258+2}, {301-2, 301+2},
                    {344-2, 344+2}}},
    {"human60-55", {{33-2, 33+2}, {66-2, 66+2}, {99-2, 99+2}, {132-2, 132+2},
                    {165-2, 165+2}, {198-2, 198+2}, {231-2, 231+2},
                    {264-2, 264+2}}},
    {"human60-75", {{23-2, 23+2}, {46-2, 46+2}, {69-2, 69+2}, {92-2, 92+2},
                    {115-2, 115+2}, {138-2, 138+2}, {161-2, 161+2},
                    {184-2, 184+2}}}
};


int main(int argc, char *argv[]) {
    BWT assembler(argv[1], std::stoi(argv[2]), std::stoi(argv[3]), argv[4]);
    if (argc >= 6) {
        std::string input(argv[5]);
        if (coverage_regions.find(input) != coverage_regions.end()) {
            std::vector<std::vector<int>> si_limits = coverage_regions[input];
            if (argc == 7) {
                for (int i = si_limits.size()-1; i >= 0; i--) {
                    std::cout << "\n----------------\n"
                              << "Starting coverage region "
                              << si_limits[i][0] << " to "
                              << si_limits[i][1] << "\n";
                    assembler.set_max_min_coverage(si_limits[i][0],
                                                   si_limits[i][1]);
                    if (input == "human30-35" && i == 6) {
                        assembler.generate_all_contigs(8);
                    } else if (input == "human30-15" && i == 6) {
                        assembler.generate_all_contigs(8);
                    } else {
                        assembler.generate_all_contigs(i+1);
                    }
                }
            } else {
                for (int i = 0; i < si_limits.size(); i++) {
                    std::cout << "\n----------------\n"
                              << "Starting coverage region "
                              << si_limits[i][0] << " to "
                              << si_limits[i][1] << "\n";
                    assembler.set_max_min_coverage(si_limits[i][0],
                                                   si_limits[i][1]);
                    if (input == "human30-35" && i == 6) {
                        assembler.generate_all_contigs(8);
                    } else if (input == "human30-15" && i == 6) {
                        assembler.generate_all_contigs(8);
                    } else {
                        assembler.generate_all_contigs(i+1);
                    }                }
            }
        } else {
            std::cout << "Wrong input" << std::endl;
            exit(1);
        }
    } else { assembler.generate_all_contigs(); }

    std::cout << std::endl;
}
