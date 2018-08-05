// Copyright [2018] Isla Carson

#define BOOST_TEST_MODULE AssemblerTests
#include <cmath>
#include <time.h>
#include <fstream>
#include <unordered_map>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include "BWT.hpp"

void is_reverse_complement(const std::string &forward,
                           const std::string &reverse);

BOOST_AUTO_TEST_SUITE(BurrowsWheelerTransform_and_FMIndex);

    BOOST_AUTO_TEST_SUITE(Initialisation);


        BOOST_AUTO_TEST_CASE(TotalCounts, *utf::enabled()) {
            /*
            *  Ensure cumulative count of the BWT characters matches the last row.
            */

            BWT test_bwt("../testdata/BWT/Enchanced_BWT_1_read.txt");
            BOOST_CHECK_MESSAGE(test_bwt.get_length() + 1 == 228,
                                "The Enchanced SA length is "
                                << test_bwt.get_length() << "+1, not 228.");

            BOOST_CHECK_MESSAGE(test_bwt.get_total_smaller_than_class(0) == 0,
                                "Total smaller than '.' should be 0.");
            // Cumulatively total last row of m_cumulative_count
            int last_row = test_bwt.get_length();
            int running_total = 0;
            for (int i = 0; i < 5; i++) {
                running_total += test_bwt.get_cumulative_count(last_row, i);
                BOOST_CHECK_MESSAGE(running_total ==
                                    test_bwt.get_total_smaller_than_class(i+1),
                                "running total of cumulative counts is "
                                << running_total
                                << " and total_smaller_than is "
                                << test_bwt.get_total_smaller_than_class(i+1));
            }
        }



        BOOST_AUTO_TEST_CASE(NumberOfTerminationChars, *utf::enabled()) {
            /*
            *  Check that the number of words in the dictionary is the same as the
            *  number of termination characters '`' in the BWT.
            *  Note that all termination chars are represented at the start of the BWT.
            */
            BWT test_bwt("../testdata/BWT/Enchanced_BWT_1_read.txt");

            Suffix_Interval test_si = test_bwt.new_SI();

            BOOST_CHECK_MESSAGE(test_si.lower == 0,
                                "The lower SI bound is "
                                << test_si.lower << ", not 0.");
            BOOST_CHECK_MESSAGE(test_si.upper == 226,
                                "The upper SI bound, " << test_si.upper
                                << ", is not the number of words in the "
                                << "dictionary, 227 (-1) as 0 indexed.");

            for (int i = 0; i < 5; i++) {
                Suffix_Interval count = test_bwt.update_backwards(test_si, i);
                BOOST_CHECK_MESSAGE(count.lower ==
                                    test_bwt.get_total_smaller_than_class(i),
                                 "The lower SI bound is "
                                << count.lower << ", not "
                                << test_bwt.get_total_smaller_than_class(i)
                                << ".");

                BOOST_CHECK_MESSAGE(count.size() ==
                    test_bwt.get_cumulative_count(test_bwt.get_length(), i),
                    "The size of the SI is " << count.size() << " not "
                    << test_bwt.get_cumulative_count(test_bwt.get_length(), i)
                    << ", which is the number of chars with eq class " << i);
            }
        }


    BOOST_AUTO_TEST_SUITE_END();

    BOOST_AUTO_TEST_SUITE(StringMatching);

        BOOST_AUTO_TEST_CASE(MersInTheDictAreInTheDict, *utf::enabled()) {
            /*
            *  Check that k-mers in the dictionary, are in the SA.
            */

            const char* path = "../testdata/BWT/Enchanced_BWT_1_read.txt";
            BWT test_bwt(path);
            std::vector<std::string> frags_in_dict = {
                "agtTAGGactattcgaacattatgtcacaaacgtgatgtcacaaagccgaattgtctggagttaagactatacGAACAttatgaaacaaacgtgatgtcac",
                "attgggcaCAGACGGAGTAGGGCAGCCTTACGTACAGATACAGATACAAACGAGAGACCAAATCATAACAGCCAACAATttggcacagacggagtagggca",
                "attgg",
                "aaaaa",
                "cCcc",
                "tttTt"};

            std::vector<std::string> suffices_in_dict;
            for (auto mer : frags_in_dict) {
                for (int i = 0; i < mer.length() - 2; i++)
                suffices_in_dict.push_back(mer.substr(i));
            }

            for (auto mer : suffices_in_dict) {
                Suffix_Interval current_si = test_bwt.new_SI();
                int i = mer.length() - 1;

                // Backwards search through the mer, terminating if not found.
                while (current_si.lower <= current_si.upper && i >= 0) {
                    int char_eq = test_bwt.convert_base_to_class(mer[i]);
                    current_si = test_bwt.update_backwards(current_si, char_eq);
                    i -= 1;
                }

                BOOST_CHECK_MESSAGE(current_si.lower <= current_si.upper,
                    "\n" << mer << " is in dictionary, but not found during "
                    << "manual search with an SI of "
                    << current_si.lower << ", " << current_si.upper);

                BOOST_CHECK_MESSAGE(test_bwt.contains_read(mer),
                    "\n" << mer << " is in dictionary, but not found during "
                    << "automatic search with an SI of "
                    << current_si.lower << ", " << current_si.upper);
            }
        }


        BOOST_AUTO_TEST_CASE(AbsentMersAreNOTInTheDict, *utf::enabled()) {
            /*
            *  Check that k-mers in the dictionary, are in the SA.
            */

            const char* path = "../testdata/BWT/Enchanced_BWT_1_read.txt";
            BWT test_bwt(path);
            std::vector<std::string> frag_not_in_dict = {
                "attgag",
                "aaaaaa",
                "cCccc",
                "tttTtt",
                "acgtacgt",
                "tcgtcgat"};

            for (auto mer : frag_not_in_dict) {
                Suffix_Interval current_si = test_bwt.new_SI();
                int i = mer.length() - 1;

                // Backwards search through the mer, terminating if not found.
                while (current_si.lower <= current_si.upper && i >= 0) {
                    int char_eq = test_bwt.convert_base_to_class(mer[i]);
                    current_si = test_bwt.update_backwards(current_si, mer[i]);
                    i -= 1;
                }

                BOOST_CHECK_MESSAGE(current_si.lower > current_si.upper,
                    "\n" << mer << " is not in the reads, but was found "
                    << "during manual search with an SI of "
                    << current_si.lower << ", " << current_si.upper);

                BOOST_CHECK_MESSAGE(!test_bwt.contains_read(mer),
                    "\n" << mer << " is not in the reads, but was found"
                    << " during automatic search");
            }
        }

        BOOST_AUTO_TEST_CASE(ReverseSISizing, *utf::enabled()) {
            const char* path = "../testdata/BWT/10000_reads__min_frag_4__sort_depth_8_with_rc.fm";
            BWT test_bwt(path);
            std::vector<std::string> eight_mers = {"ACATTATG", "ACAGATAC",
            "TTCCATTC", "AGTTTGAG", "TTGCTAGC", "AAAGTAAG", "CAGTACTG",
            "TTAGGACT", "CCCCCAAC", "AGCTGTTT", "AAAGAGCT", "TCAATGCT",
            "AAACGAAG", "CGTTTGCA", "AGAGACAA", "TGTCTATG", "GGAGGAAG",
            "CATATCTG", "AAAGAGTG", "ATCGGTCG"};

            for (auto mer : eight_mers) {
                Suffix_Interval normal_si = test_bwt.new_SI();
                int i = mer.length() - 1;

                // For each, perform a backwards search and a search for the
                //  rc.  Check the SI sizes are equal.
                while (normal_si.lower <= normal_si.upper && i >= 0) {
                    normal_si = test_bwt.update_backwards(normal_si, mer[i]);
                    i -= 1;
                }

                std::string reverse_rc = reverse_complement(mer);
                Suffix_Interval rc_si = test_bwt.new_SI();
                i = reverse_rc.length() - 1;

                // For each, perform a backwards search and a search for the
                //  rc.  Check the SI sizes are equal.
                while (rc_si.lower <= rc_si.upper && i >= 0) {
                    rc_si = test_bwt.update_backwards(rc_si, reverse_rc[i]);
                    i -= 1;
                }

                BOOST_CHECK_MESSAGE(normal_si.valid() == rc_si.valid(),
                                    "One Si is valid and the other is not.");
                if (normal_si.valid()) {
                    BOOST_CHECK_MESSAGE(normal_si.size() == rc_si.size(),
                        "normal and reverse si sizes do not match.");
                }
            }
        }



    BOOST_AUTO_TEST_SUITE_END();

    BOOST_AUTO_TEST_SUITE(LCPFlags);

        BOOST_AUTO_TEST_CASE(LCPFlagBoundaries, *utf::enabled()) {
            const char* path = "../testdata/BWT/10000_reads__min_frag_4__sort_depth_8.txt";
            BWT test_bwt(path);
            std::vector<std::string> eight_mers = {"ACATTATG", "ACAGATAC",
            "TTCCATTC", "AGTTTGAG", "TTGCTAGC", "AAAGTAAG", "CAGTACTG",
            "TTAGGACT", "CCCCCAAC", "AGCTGTTT", "AAAGAGCT", "TCAATGCT",
            "AAACGAAG", "CGTTTGCA", "AGAGACAA", "TGTCTATG", "GGAGGAAG",
            "CATATCTG", "AAAGAGTG", "ATCGGTCG"};

            for (auto mer : eight_mers) {
                Suffix_Interval current_si = test_bwt.new_SI();
                int i = mer.length() - 1;

                // For each, perform a backwards search and assert that each SI
                // does not intersect any LCP classes (i.e. it contains >=1
                // full LCP classes.)
                while (current_si.lower <= current_si.upper && i >= 0) {
                    // Lower boundary is an LCP class boundary
                    if (current_si.lower > 0) {
                        BOOST_CHECK_MESSAGE(
                            test_bwt.get_lcp_flag(current_si.lower-1)
                                != test_bwt.get_lcp_flag(current_si.lower),
                            "Lower boundary crosses an lcp class.");
                    }
                    // Upper boundary is an LCP class boundary
                    BOOST_CHECK_MESSAGE(test_bwt.get_lcp_flag(current_si.upper)
                                != test_bwt.get_lcp_flag(current_si.upper+1),
                            "Upper boundary crosses an lcp class.");

                    current_si = test_bwt.update_backwards(current_si, mer[i]);
                    i -= 1;
                }

                // The SI at length 8 (if it exists) is a single LCP class.
                if (current_si.lower <= current_si.upper) {
                    // Lower boundary is an LCP class boundary
                    if (current_si.lower > 0) {
                        BOOST_CHECK_MESSAGE(
                            test_bwt.get_lcp_flag(current_si.lower-1)
                                != test_bwt.get_lcp_flag(current_si.lower),
                            "Lower boundary crosses an lcp class.");
                    }

                    // Upper boundary is an LCP class boundary
                    BOOST_CHECK_MESSAGE(test_bwt.get_lcp_flag(current_si.upper)
                                != test_bwt.get_lcp_flag(current_si.upper+1),
                                "Upper boundary crosses an lcp class.");

                    // No class boundaries contained within.
                    for (int i = current_si.lower; i < current_si.upper; i++) {
                        BOOST_CHECK_MESSAGE(test_bwt.get_lcp_flag(i)
                                                == test_bwt.get_lcp_flag(i+1),
                            "Class contains more than one LCP class detected"
                            << " at " << i);
                    }
                }
            }
        }

        BOOST_AUTO_TEST_CASE(LCPFlagLetterDrop, *utf::enabled()) {
            const char* path = "../testdata/BWT/10000_reads__min_frag_4__sort_depth_8.txt";
            BWT test_bwt(path);

            // Create a bunch of 9-mers
            std::vector<std::string> nine_mers = {"ATGTGTGAT", "GAAAAGTAA",
                "GACGTGTGT", "TTGTCTTGA", "ATATCACCT", "CTCCGTCTG",
                "AACCAAATC", "CCTTGATAA", "GCTCACCAA", "GCACAGACG",
                "TGTGCAAAT", "GGGCAGGGG", "CAGGAGGTG", "TTGGATCTG",
                "CTATTAGGA", "CACGCACGC", "AACACACTA", "TCTCTGCTT",
                "GACCAAATC", "TGTGTTGTT", "ACCAGCATT"};

            for (auto mer : nine_mers) {
                // Find the SI for the first 8 bases, and then
                // seperately for the full 9-mer.
                Suffix_Interval prefix_si = test_bwt.new_SI();
                int i = 7;  // To exclude the last letter.
                while (prefix_si.lower <= prefix_si.upper && i >= 0) {
                    prefix_si = test_bwt.update_backwards(prefix_si, mer[i]);
                    i -= 1;
                }

                Suffix_Interval all_si = test_bwt.new_SI();
                i = 8;  // To exclude the last letter.
                while (all_si.lower <= all_si.upper && i >= 0) {
                    all_si = test_bwt.update_backwards(all_si, mer[i]);
                    i -= 1;
                }

                // Use the LCP flags to scan outwards and find the SI for the
                // 9-mer with the last letter dropped off.
                all_si = test_bwt.drop_last_base(all_si);

                // Test that this is the same SI as that found for the 8-mer
                // prefix
                BOOST_CHECK_MESSAGE(prefix_si.lower == all_si.lower,
                    "Prefix SI has lower bound " << prefix_si.lower
                    << " whilst all_SI is " << all_si.lower);
                BOOST_CHECK_MESSAGE(prefix_si.upper == all_si.upper,
                    "Prefix SI has upper bound " << prefix_si.upper
                    << " whilst all_SI is " << all_si.upper);
                BOOST_CHECK_MESSAGE(prefix_si.match_length
                                                        == all_si.match_length,
                    "Prefix SI has match_length " << prefix_si.match_length
                    << " whilst all_SI has " << all_si.match_length);
            }

            Suffix_Interval hits_0 = test_bwt.drop_last_base({0, 9977, 0});
            BOOST_CHECK_MESSAGE(hits_0.lower == 0,
                    "hits_0 SI has lower bound " << hits_0.lower
                    << " but should be " << 0);
                BOOST_CHECK_MESSAGE(hits_0.upper == 9999,
                    "hits_0 SI has upper bound " << hits_0.upper
                    << " but should be " << 9999);

            Suffix_Interval first_cls = test_bwt.drop_last_base({31, 9900, 0});
            BOOST_CHECK_MESSAGE(first_cls.lower == 0,
                    "first_cls SI has lower bound " << first_cls.lower
                    << " but should be " << 0);
                BOOST_CHECK_MESSAGE(first_cls.upper == 9999,
                    "first_cls SI has upper bound " << first_cls.upper
                    << " but should be " << 9999);

            Suffix_Interval full_cls = test_bwt.drop_last_base({10000, 12356,
                                                                0});
            BOOST_CHECK_MESSAGE(full_cls.lower == 10000,
                    "full_cls SI has lower bound " << full_cls.lower
                    << " but should be " << 10000);
                BOOST_CHECK_MESSAGE(full_cls.upper == 12356,
                    "full_cls SI has upper bound " << full_cls.upper
                    << " but should be " << 12356);

            Suffix_Interval last_cls = test_bwt.drop_last_base({956289, 956374,
                                                                0});
            BOOST_CHECK_MESSAGE(last_cls.lower == 956289,
                    "last_cls SI has lower bound " << last_cls.lower
                    << " but should be " << 956289);
                BOOST_CHECK_MESSAGE(last_cls.upper == 956374,
                    "last_cls SI has upper bound " << last_cls.upper
                    << " but should be " << 956374);
        }


    BOOST_AUTO_TEST_SUITE_END();

    BOOST_AUTO_TEST_SUITE(MatchBehaviour);

        BOOST_AUTO_TEST_CASE(ContigsComplement, *utf::enabled()) {

            const char* path = "../testdata/BWT/10000_reads__min_frag_4__sort_depth_8_with_rc.fm";
            BWT test_bwt(path, 16, 400);

            // Create a bunch of 9-mers
            std::vector<std::string> eight_mers = {
                "ATGTGTGA", "GAAAAGTA",
                "GACGTGTG", "TTGTCTTG", "ATATCACC", "CTCCGTCT",
                "AACCAAAT", "CCTTGATA", "GCTCACCA", "GCACAGAC",
                "TGTGCAAA", "GGGCAGGG", "CAGGAGGT", "TTGGATCT",
                "CTATTAGG"
                , "CACGCACG", "AACACACT", "TCTCTGCT",
                "GACCAAAT", "TGTGTTGT", "ACCAGCAT"
            };

            for (auto mer : eight_mers) {
                std::string contig = test_bwt.find_contig(mer);
                test_bwt.clear_visited_log();
                std::string rc = reverse_complement(mer);
                std::string revcon = test_bwt.find_contig(rc);

                is_reverse_complement(contig, revcon);
            }
        }

        BOOST_AUTO_TEST_CASE(ComplementVisited, *utf::enabled()) {

            const char* path = "../testdata/BWT/10000_reads__min_frag_4__sort_depth_8_with_rc.fm";
            BWT test_bwt(path, 40, 360);

            // Create a bunch of 9-mers
            std::vector<std::string> eight_mers = {
                "AAAACTAG",
                "AAAACTGA",
                "AAAATAAA",
                "AAAATATG",
                "AAAGTTGT",
                "AACACTAT",
                "AACCAAAT",
                "AACACACT",
                "CCTTGATA",
                "AAAAAACA",  //  aAAAAAACAa
                "ATGAACAA"  //  ATGAACAAag
            };

            for (auto mer : eight_mers) {
                test_bwt.clear_visited_log();
                std::string contig = test_bwt.find_contig(mer);
                std::string rc_contig = reverse_complement(contig);
                for (int i = 0;
                    i < contig.length()-test_bwt.get_sort_depth();
                    i++) {
                    Suffix_Interval current_si = test_bwt.get_SI_for(
                                contig.substr(i, test_bwt.get_sort_depth()));

                    bool is_visited = test_bwt.get_visited(current_si.lower);
                    bool bad_size =
                            (current_si.size() < test_bwt.get_min_coverage())
                         || (current_si.size() > test_bwt.get_max_coverage());

                    //  Valid states:
                        // Visited and good size.
                        // Not visited and bad size.

                    // Incorrect states
                        // Visited & bad size.
                        // Not visited and good size

                    std::string error_msg;
                    if (is_visited) {
                        error_msg = "SI has been visited and has bad size.";
                    } else {
                        error_msg = "SI has not been visited and has good size.";
                    }

                    BOOST_CHECK_MESSAGE((is_visited && !bad_size)
                                        || (!is_visited && bad_size),
                        "BWT should be visited at " << current_si.lower
                        << " with kmer " << mer
                        << ", contig " << contig << " and contig substr "
                        << contig.substr(i, test_bwt.get_sort_depth())
                        << ", but isn't.");

                    current_si = test_bwt.get_SI_for(
                            rc_contig.substr(i, test_bwt.get_sort_depth()));
                    is_visited = test_bwt.get_visited(current_si.lower);
                    bad_size =
                            (current_si.size() < test_bwt.get_min_coverage())
                         || (current_si.size() > test_bwt.get_max_coverage());

                    BOOST_CHECK_MESSAGE(!test_bwt.is_valid(current_si),
                        "BWT should be visited at " << current_si.lower
                        << " with kmer " << mer
                        << ", RCcontig " << rc_contig
                        << " and RCcontig substr "
                        << rc_contig.substr(i, test_bwt.get_sort_depth())
                        << ", but isn't.");
                }
            }
        }

        BOOST_AUTO_TEST_CASE(AllSIsVisitedOrWrongSize, *utf::enabled()) {
            const char* path = "../testdata/BWT/10000_reads__min_frag_4__sort_depth_16_with_rc.fm";
            BWT test_bwt(path, 20, 35, "../testdata/results/");

            // {45, 55}, {90, 110}

            test_bwt.generate_all_contigs();
            std::bitset<160> kmer_seed(0);
            uint_least64_t counter = 0;
            while (test_bwt.max_seed_passed(kmer_seed)) {
                Suffix_Interval forward = test_bwt.get_SI_for_seed(kmer_seed);

                if (forward.valid()) {
                    //  Valid states:
                        // Visited and good size.
                        // Not visited and bad size.

                    // Incorrect states
                        // Visited & bad size.
                        // Not visited and good size

                    bool is_visited = test_bwt.get_visited(forward.lower);
                    bool bad_size =
                                (forward.size() < test_bwt.get_min_coverage())
                            || (forward.size() > test_bwt.get_max_coverage());

                    std::string error_msg;
                    if (is_visited) {
                        error_msg = "SI has been visited and has bad size.";
                    } else {
                        error_msg = "SI hasn't been visited and has good size";
                    }

                    BOOST_CHECK_MESSAGE((is_visited && !bad_size)
                                            || (!is_visited && bad_size),
                                        error_msg << " " << is_visited
                                        << " " << bad_size
                                        << " "
                                        <<  test_bwt.convert_to_kmer(kmer_seed)
                                        << " " << forward.size());
                }

                counter++;
                if (counter%1000000 == 0) { std::cout << counter/1000000
                                                      << "t " << std::flush; }

                // std::cout << test_bwt.convert_to_kmer(kmer_seed) << " "
                          // << forward.match_length <<"\n";
                test_bwt.add_one_at_index(0, &kmer_seed);
            }
            std::cout << "Final seed:" << test_bwt.convert_to_kmer(kmer_seed)
                      << "  " << kmer_seed.to_string();

        }


    BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();



void is_reverse_complement(const std::string &forward,
                           const std::string &reverse) {
    int length = forward.length();
    BOOST_CHECK_MESSAGE(length == reverse.length(),
                "Forward and reverse complement contig lengths do not match.");

    bool match = true;
    for (int i = 0; i < length; i++) {
        if (forward[i] != complement(reverse[length-i-1])) {
            match = false;
        }
    }
    BOOST_CHECK_MESSAGE(match,
                        "Forward ("<< forward
                        << ") and reverse do not match, with forward:\n"
                        << forward
                        << "\n" <<
                        reverse_complement(reverse) << "\n");
}
