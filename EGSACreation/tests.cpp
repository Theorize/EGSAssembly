// Copyright [2018] Isla Carson

#define BOOST_TEST_MODULE ConstructSuffixArrayTests
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;


#include "ConstructSA.hpp"
#include "LinearSA.hpp"

template <class SA> void test_order(const SA &test_SA,
                                       const std::string &concat_reads);
int path_inverse(std::stack<bool> path);
std::string get_cyclic_shift_at(int i, int length,
                                const std::string &sequence);
template <class SA> void try_full_alg(std::string file_name,
                        std::unordered_map<std::string, uint_least64_t> input,
                        SA &test_SA);
std::string dump_reads_to_string(std::string file_name, int file_length);
void compare_prob_desc(ProblemDescription first, ProblemDescription second);

BOOST_AUTO_TEST_SUITE(ProbDesc);

BOOST_AUTO_TEST_CASE(TwoReadsStats) {
    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 5;
    input["-min_frag_len"] = 5;
    input["-precise_depth"] = 1;
    ProblemDescription stats_only("../testdata/reads-1.fastq", input);

    BOOST_CHECK_MESSAGE(stats_only.get_SA_length() == 222,
                          "Problem length is "
                          << stats_only.get_SA_length() << ", not 222");
    BOOST_CHECK_MESSAGE(stats_only.get_depth() == 5,
                        "get_depth fn is "<< stats_only.get_depth()
                        << ", not 5");
    BOOST_CHECK_MESSAGE(stats_only.get_min_frag_len() == 5,
                        "get_min_frag_len fn is "
                        << stats_only.get_min_frag_len()
                        << ", not 5");
    BOOST_CHECK_MESSAGE(stats_only.get_precision_req() == 1,
                        "get_precision_req fn is "
                        << stats_only.get_precision_req()
                        << ", not 1");
    BOOST_CHECK_MESSAGE(stats_only.get_read_number() == 3,
                        "get_read_number fn is "
                        << stats_only.get_read_number() << ", not 3");
    BOOST_CHECK_MESSAGE(stats_only.get_percentage_useful_reads() == 100.0,
                        "get_percentage_useful_reads fn is "
                        << stats_only.get_percentage_useful_reads()
                        << ", not 100.0");
    BOOST_CHECK_MESSAGE(stats_only.get_fragment_number() == 5,
                        "get_fragment_number fn is "
                        << stats_only.get_fragment_number()
                        << ", not 5");
    BOOST_CHECK_MESSAGE(stats_only.get_avg_frag_length() == 43.4,
                        "get_avg_frag_length fn is "
                        << stats_only.get_avg_frag_length()
                        << ", not 43.4");
    BOOST_CHECK_MESSAGE(stats_only.get_total_bases() == 224,
                        "get_total_bases fn is "
                        << stats_only.get_total_bases()
                        << ", not 224");
    BOOST_CHECK_MESSAGE(stats_only.get_num_useable_bases() == 217,
                        "get_num_useable_bases fn is "
                        << stats_only.get_num_useable_bases()
                        << ", not 217");
    BOOST_CHECK_MESSAGE(stats_only.get_percentage_useful_bases() == 96.875,
                        "get_percentage_useful_bases fn is "
                        << stats_only.get_percentage_useful_bases()
                        << ", not 96.875");
}

BOOST_AUTO_TEST_CASE(ProbDescFnOutputsMatch) {
    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 5;
    input["-min_frag_len"] = 5;

    std::vector<std::string> paths {"../testdata/reads-1.fastq",
                                   "../testdata/reads-100ish.fastq",
                                   "../testdata/reads-10000.fastq"};

    for (auto file : paths) {
        ProblemDescription stats_only(file, input);

        ProblemDescription to_file(file,
                                   "../testdata/Reads.txt", input);
        std::string reads = dump_reads_to_string("../testdata/Reads.txt",
                                                 to_file.get_SA_length());

        ProblemDescription to_classes(file, input, 1);
        std::vector<uint_least64_t> classes =
                                        to_classes.process_problem_to_vector();

        compare_prob_desc(stats_only, to_file);
        compare_prob_desc(stats_only, to_classes);
        compare_prob_desc(to_file, to_classes);

        for (int i = 0; i < to_file.get_SA_length(); i++) {
            char c = '.';
            if (classes[i] == 1) { c = 'A';
            } else if (classes[i] == 2) { c = 'C';
            } else if (classes[i] == 3) { c = 'G';
            } else if (classes[i] == 4) { c = 'T'; }
            BOOST_CHECK_MESSAGE(reads[i] == c,
                                "to_file and classes do not match"
                                << " at "  << i << ". to_file: " << reads[i]
                                << " classes: " << classes[i] << " ie, "<< c);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


BOOST_AUTO_TEST_SUITE(ThreadedAlg);

BOOST_AUTO_TEST_SUITE(Light);

BOOST_AUTO_TEST_CASE(ResolvesToCorrectDepth, *utf::label("Initialise")) {

    for (int d = 5; d <= 100; d++) {
        std::unordered_map<std::string, uint_least64_t> input;
        input["-sort_depth"] = d;
        input["-precise_depth"] = 1;
        ConstructSA test_problem("../testdata/reads-100ish.fastq",
                                        input);
        std::stack<bool> path = test_problem.copy_alg_steps();

        BOOST_CHECK_MESSAGE(path.top() == 1, "Path top not 1 for d = " << d);

        BOOST_CHECK_MESSAGE(d == path_inverse(path),
                            "Path does not result in required depth = " << d);
    }
     // std::cout << " ... Done\n\n" << std::flush;
}

BOOST_AUTO_TEST_CASE(CorrectCharIndices, *utf::label("Initialise")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;
    // Create test class
    ProblemDescription test_PD("../testdata/reads-100ish.fastq",
                               "../testdata/Reads.txt",
                               input);
    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    // Count the number of 'T's in the reads.txt file
    std::ifstream is("../testdata/Reads.txt");

    char c;
    int t_count = 0;
    while (is.get(c)) {
        // Get the order of the char, and increase by 1.
        if (c == 'T') { t_count += 1; }
    }
    is.close();
    // Verify the number of bases in the file is the same as the start index
    // of 'T' plus the number of 'T's.
    int total = t_count + test_SA.get_base_index(4);
    BOOST_CHECK_MESSAGE(total == test_SA.get_SA_length(),
        "Number of 'T's: " << t_count
        << "\n Number of bases smaller than T: " << test_SA.get_base_index(4));
    // std::cout << " ... Done\n\n" << std::flush;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


BOOST_AUTO_TEST_CASE(CorrectMinClassForBases, *utf::label("SortSetup")) {
        std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;

    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    test_SA.initialise_alg();
    // For each base.
    for (int run = 0; run < 4; run++) {
        std::cout.setstate(std::ios_base::failbit);
        test_SA.sort_suffix_pairs(false);
        test_SA.update_classes(0);
        std::cout.clear();
        for (int i=0; i < 5; i++) {
            BOOST_CHECK_MESSAGE(test_SA.get_min_class_for_base(i) ==
                      test_SA.get_class_at(test_SA.get_suffix_at(
                                                test_SA.get_base_index(i))),
                      "i is: " << i << " Iteration Number is: " << run);
        }
    }
    // std::cout << " ... Done\n\n" << std::flush;
}


BOOST_AUTO_TEST_CASE(CyclicShiftsLen1FromFile, *utf::label("SortSetup")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;

    ProblemDescription test_PD_with_file("../testdata/reads-100ish.fastq",
                                         "../testdata/Reads.txt",  input);
    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length() == test_SA.get_SA_length());

    test_SA.initialise_alg();

    char bases[5] = {'.', 'A', 'C', 'G', 'T'};
    std::map<char, int[2]> suffix_interval;
    suffix_interval['.'][0] = 0;
    suffix_interval['.'][1] = test_SA.get_base_index(0);
    for (int i=1; i < 4; i++) {
        suffix_interval[bases[i]][0] = test_SA.get_base_index(i-1);
        suffix_interval[bases[i]][0] = test_SA.get_base_index(i);
    }
    suffix_interval['T'][0] = test_SA.get_base_index(4);
    suffix_interval['T'][1] = test_SA.get_SA_length();

    for (auto c : bases) {
        // Check that each base is in the correct position.
        for (int i = suffix_interval[c][0]; i < suffix_interval[c][1]; i++) {
            BOOST_CHECK_MESSAGE(concat_reads[test_SA.get_suffix_at(i)] == c,
                "position in the SA: " << i
                << "\nvalue in the SA: " <<  test_SA.get_suffix_at(i));
        }
    }
    // std::cout << " ... Done\n\n" << std::flush;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


BOOST_AUTO_TEST_CASE(ComparisonDoubleSortDepth, *utf::label("CompSort")) {
        std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;

    ProblemDescription test_PD_with_file("../testdata/reads-100ish.fastq",
                                         "../testdata/Reads.txt", input);
    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length() == test_SA.get_SA_length());

    test_SA.initialise_alg();

    for (int j = 0; j < 5; j++) {
        // Testing Doubling Functionality.  So mode = 0.
        bool mode = 0;
        bool choice = 0;
        std::cout.setstate(std::ios_base::failbit);
        test_SA.sort_suffix_pairs(mode);
        test_SA.update_classes(mode);
        test_SA.set_all_sort_choices(choice);
        std::cout.clear();
        test_order(test_SA, concat_reads);
    }
    // std::cout << " ... Done\n\n" << std::flush;
}


BOOST_AUTO_TEST_CASE(ComparisonAdd1ToSortDepth, *utf::label("CompSort")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;
    ProblemDescription test_PD_with_file("../testdata/reads-100ish.fastq",
                                         "../testdata/Reads.txt", input);
    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length() == test_SA.get_SA_length());
    test_SA.initialise_alg();

    for (int j = 0; j < 5; j++) {
        // Testing Add 1 Functionality.  So mode = 1.
        bool mode = 1;
        bool choice = 0;
        std::cout.setstate(std::ios_base::failbit);
        test_SA.sort_suffix_pairs(mode);
        test_SA.update_classes(mode);
        test_SA.set_all_sort_choices(choice);
        std::cout.clear();
        test_order(test_SA, concat_reads);
    }
    // std::cout << " ... Done\n\n" << std::flush;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


BOOST_AUTO_TEST_CASE(CombDoubleSortDepth, *utf::label("CombinationSort")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;

    ProblemDescription test_PD_with_file("../testdata/reads-100ish.fastq",
                                         "../testdata/Reads.txt",  input);
    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length() == test_SA.get_SA_length());

    test_SA.initialise_alg();

    for (int j = 0; j < 5; j++) {
        // Testing Doubling Functionality.  So mode = 0.
        bool mode = 0;
        bool choice = 1;
        std::cout.setstate(std::ios_base::failbit);
        test_SA.sort_suffix_pairs(mode);
        test_SA.update_classes(mode);
        test_SA.set_all_sort_choices(choice);
        std::cout.clear();
        test_order(test_SA, concat_reads);
    }
    // std::cout << " ... Done\n\n" << std::flush;
}


BOOST_AUTO_TEST_CASE(CombAdd1ToSortDepth, *utf::label("CombinationSort")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;

    ProblemDescription test_PD_with_file("../testdata/reads-100ish.fastq",
                                         "../testdata/Reads.txt", input);
    ConstructSA test_SA("../testdata/reads-100ish.fastq", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length(), test_SA.get_SA_length());
    test_SA.initialise_alg();

    for (int j = 0; j < 5; j++) {
        // Testing Add 1 Functionality.  So mode = 1.
        bool mode = 1;
        bool choice = 1;
        std::cout.setstate(std::ios_base::failbit);
        test_SA.sort_suffix_pairs(mode);
        test_SA.update_classes(mode);
        test_SA.set_all_sort_choices(choice);
        std::cout.clear();
        test_order(test_SA, concat_reads);
    }
    // std::cout << " ... Done\n\n" << std::flush;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(FullAlgRun, *utf::label("FullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    input["-sort_depth"] = 20;
    input["-precise_depth"] = 0;
    input["-alg_limit"] = 800;
    input["-min_frag_len"] = 75;
    input["-num_threads"] = 2;
    std::string file_name("../testdata/reads-100ish.fastq");
    ConstructSA test_SA(file_name, input);

    try_full_alg(file_name, input, test_SA);

    // std::cout << " ... Done\n\n" << std::flush;
}


BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(VariedParams);

BOOST_AUTO_TEST_CASE(MinValues, *utf::label("FullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    // Initialise all to minimum allowed values
    input["-sort_depth"] = 1;
    input["-precise_depth"] = 1;
    input["-alg_limit"] = 1;
    input["-min_frag_len"] = 1;
    input["-num_threads"] = 1;
    std::string file_name("../testdata/reads-100ish.fastq");

    // Test
    ConstructSA precise_test_SA(file_name, input);
    try_full_alg(file_name, input, precise_test_SA);

    // Precision not required
    input["-precise_depth"] = 0;

    // Test
    ConstructSA test_SA(file_name, input);
    try_full_alg(file_name, input, test_SA);

    // std::cout << " ... Done\n\n" << std::flush;
}

BOOST_AUTO_TEST_CASE(DepthAndPrecisionVariation, *utf::label("FullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    input["-alg_limit"] = 8000;
    input["-min_frag_len"] = 20;
    input["-num_threads"] = 4;
    std::string file_name("../testdata/reads-100ish.fastq");
    for (int depth = 1; depth < 100; depth+=10) {
        // Initialise
        input["-sort_depth"] = depth;
        input["-precise_depth"] = 1;
        // Test
        ConstructSA precise_test_SA(file_name, input);
        try_full_alg(file_name, input, precise_test_SA);

        // Initialise
        input["-precise_depth"] = 0;
        // Test
        ConstructSA test_SA(file_name, input);
        try_full_alg(file_name, input, test_SA);
    }
    // std::cout << " ... Done\n\n" << std::flush;
}

BOOST_AUTO_TEST_CASE(AlgLimitVariation, *utf::label("FullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    input["-min_frag_len"] = 20;
    input["-num_threads"] = 4;
    input["-sort_depth"] = 16;
    std::string file_name("../testdata/reads-10000.fastq");
    std::cout << "    Alg limit at: "  << std::flush;
    for (int alg_limit = 1; alg_limit < 80000; alg_limit*=2) {
        std::cout << alg_limit  << ", "  << std::flush;
        // Initialise
        input["-alg_limit"] = alg_limit;
        input["-precise_depth"] = 1;
        // Test
        ConstructSA precise_test_SA(file_name, input);
        try_full_alg(file_name, input, precise_test_SA);

        // Initialise
        input["-precise_depth"] = 0;
        // Test
        ConstructSA test_SA(file_name, input);
        try_full_alg(file_name, input, test_SA);

        // If alg_limit small, skip some iterations.
        if (alg_limit == 4) { alg_limit = 64; }
        // If alg_limit big, skip some iterations.
        if (alg_limit >= 128) { alg_limit *= 2; }
    }
    // std::cout << " ... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_CASE(MinFragLenVariation, *utf::label("FullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    input["-num_threads"] = 4;
    input["-sort_depth"] = 16;
    std::string file_name("../testdata/reads-10000.fastq");
    std::cout << "    Min Frag Len at: "  << std::flush;
    for (int min_frag_len = 1; min_frag_len < 100; min_frag_len+=20) {
        std::cout << min_frag_len  << ", "  << std::flush;
        // Initialise
        input["-min_frag_len"] = min_frag_len;
        input["-precise_depth"] = 1;
        // Test
        ConstructSA precise_test_SA(file_name, input);
        try_full_alg(file_name, input, precise_test_SA);

        // Initialise
        input["-precise_depth"] = 0;
        // Test
        ConstructSA test_SA(file_name, input);
        try_full_alg(file_name, input, test_SA);
    }
    // std::cout << "\n... Done\n\n"  << std::flush;
}


BOOST_AUTO_TEST_CASE(SmallFileSizesVSNumThreads, *utf::label("FullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 32;

    int threads[3] = {1, 4, 8};
    std::string file_name[2] = {"../testdata/reads-1.fastq",
                               "../testdata/reads-10000.fastq"};

    // For each file
    for (int i = 0; i < 2; i++) {
        std::cout << "    Fastx file: " << file_name[i]
                  << "\n        Num Threads: " << std::flush;
        //  For each number of threads.
        for (int j = 0; j < 3; j++) {
            std::cout << threads[j] << ", ";
            input["-num_threads"] = threads[j];
            ConstructSA test_SA(file_name[i], input);
            try_full_alg(file_name[i], input, test_SA);
        }
        std::cout << "\n" << std::flush;
    }
    // std::cout << "... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(Heavy);

BOOST_AUTO_TEST_CASE(LargeFileDiffNumThreads,
                     *utf::label("FullAlg")) {
        std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 32;

    int threads[3] = {1, 4, 8};
    std::string file_name = "../testdata/454/Ecoli.FLX.fna";

    //  For each number of threads.
    for (int j = 0; j < 3; j++) {
        std::cout << "    Fastx file: " << file_name
                  << " with " << threads[j] << " thread(s).\n" << std::flush;
        input["-num_threads"] = threads[j];
        ConstructSA test_SA(file_name, input);
        try_full_alg(file_name, input, test_SA);
    }
    // std::cout << "\n... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(LinearAlg);
BOOST_AUTO_TEST_SUITE(Light);

BOOST_AUTO_TEST_CASE(LinearDoubleSortDepth, *utf::label("LinearCountSort")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-sort_depth"] = 10;

    ProblemDescription test_PD_with_file("../testdata/reads-100ish.fastq",
                                         "../testdata/Reads.txt", input);
    LinearSA test_SA("../testdata/reads-100ish.fastq", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length() == test_SA.get_SA_length());


    test_SA.initialise_alg();
         std::cout << "here. " << std::flush;

    for (int j = 0; j < 5; j++) {
        test_SA.sort_suffix_pairs_count_sort();
         std::cout << "here. " << std::flush;

        test_SA.update_classes();
         std::cout << "there. " << std::flush;

        test_order(test_SA, concat_reads);
        std::cout << "everywhere. " << std::flush;

    }
    // std::cout << " ... Done\n\n"  << std::flush;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------



BOOST_AUTO_TEST_CASE(LinearRun,  *utf::label("LinearFullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    input["-sort_depth"] = 25;
    input["-min_frag_len"] = 75;
    std::string file_name("../testdata/reads-100ish.fastq");

    LinearSA test_SA(file_name, input);
    try_full_alg(file_name, input, test_SA);

    // std::cout << " ... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_CASE(LinearMinValues, *utf::label("LinearFullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    // Initialise all to minimum allowed values
    input["-sort_depth"] = 1;
    input["-min_frag_len"] = 1;
    std::string file_name("../testdata/reads-100ish.fastq");

    // Test
    LinearSA test_SA(file_name, input);
    try_full_alg(file_name, input, test_SA);

    // std::cout << " ... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_CASE(LinearDepthAndMinFragLenVariation,
                     *utf::label("LinearFullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;
    std::string file_name("../testdata/reads-100ish.fastq");

    for (int depth = 1; depth < 101; depth+=10) {
        if (depth != 1) { std::cout << "\n" << std::flush; }
        std::cout << "    Depth " << depth << " & min_frag_len: "
                  << std::flush;
        for (int min_frag_len = 10; min_frag_len < 101; min_frag_len += 10) {
            std::cout << min_frag_len << ", " << std::flush;
            // Initialise
            input["-sort_depth"] = depth;
            input["-min_frag_len"] = min_frag_len;
            // Test
            LinearSA test_SA(file_name, input);
            try_full_alg(file_name, input, test_SA);
        }
        if (depth == 1) { depth = 0; }
        if (depth == 100) {std::cout << "\n" << std::flush;}
    }
    // std::cout << "... Done\n\n"  << std::flush;
}


BOOST_AUTO_TEST_CASE(LinearDepthOnlyVariation, *utf::label("LinearFullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;

    std::string file_name("../testdata/reads-100ish.fastq");
    std::cout << "    Depth: " << std::flush;
    for (int depth = 1; depth < 101; depth+=10) {
        std::cout << depth << ", " << std::flush;
        // Initialise
        input["-sort_depth"] = depth;
        // Test
        LinearSA test_SA(file_name, input);
        try_full_alg(file_name, input, test_SA);
        if (depth == 1) { depth = 0; }
    }
    // std::cout << "... Done\n\n"  << std::flush;
}



BOOST_AUTO_TEST_CASE(LinearMinFragLenOnlyVariation,
                     *utf::label("LinearFullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;
    std::string file_name("../testdata/reads-100ish.fastq");

    std::cout << "    Min_frag_len: " << std::flush;
    for (int min_frag_len = 1; min_frag_len < 101; min_frag_len += 10) {
        std::cout << min_frag_len << ", " << std::flush;
        // Initialise
        input["-min_frag_len"] = min_frag_len;
        // Test
        LinearSA test_SA(file_name, input);
        try_full_alg(file_name, input, test_SA);
        if (min_frag_len == 1) { min_frag_len = 0; }
    }
    // std::cout << " ... Done\n\n"  << std::flush;
}


BOOST_AUTO_TEST_CASE(LinearSmallFileSizes, *utf::label("LinearFullAlg")) {

    std::unordered_map<std::string, uint_least64_t> input;
    input["-min_frag_len"] = 100;

    std::string file_name[3] = {"../testdata/reads-1.fastq",
                                "../testdata/reads-100ish.fastq",
                                "../testdata/reads-10000.fastq"};
    std::cout << "    Fastx file: " << std::flush;
    // For each file
    for (int i = 0; i < 3; i++) {
        std::cout << file_name[i] << ", " << std::flush;
        //  For each number of threads.
        LinearSA test_SA(file_name[i], input);
        try_full_alg(file_name[i], input, test_SA);
        std::cout << "\n        " << std::flush;
    }
    // std::cout << "\n... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE(Heavy);

BOOST_AUTO_TEST_CASE(LinearLargeFile, *utf::label("HeavyTests")) {

    // Setup
    std::unordered_map<std::string, uint_least64_t> input;
    input["-min_frag_len"] = 49;
    std::string file_name = "../testdata/454/Ecoli.FLX.fna";

    // Tell user
    std::cout << "    Fastx file: " << file_name << "\n" << std::flush;

    // Run test
    LinearSA test_SA(file_name, input);
    try_full_alg(file_name, input, test_SA);

    // std::cout << "\n... Done\n\n"  << std::flush;
}

BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE_END();


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


void compare_prob_desc(ProblemDescription first, ProblemDescription second) {
    BOOST_REQUIRE_MESSAGE(first.get_SA_length() ==
                                        second.get_SA_length(),
                          "Problem lengths don't match. first is "
                          << first.get_SA_length() << " second is "
                          << second.get_SA_length());
    BOOST_CHECK_MESSAGE(second.get_depth() == first.get_depth(),
                        "get_depth fn don't match."
                        << "first is " << first.get_depth()
                        << ", classes is " << second.get_depth());
    BOOST_CHECK_MESSAGE(second.get_min_frag_len()
                            == first.get_min_frag_len(),
                        "get_min_frag_len fn don't match."
                        << "first is " << first.get_min_frag_len()
                        << ", classes is " << second.get_min_frag_len());
    BOOST_CHECK_MESSAGE(second.get_precision_req()
                                == first.get_precision_req(),
                        "get_precision_req fn don't match."
                        << "first is " << first.get_precision_req()
                        << ", classes is " << second.get_precision_req());
    BOOST_CHECK_MESSAGE(second.get_read_number()
                                == first.get_read_number(),
                        "get_read_number fn don't match."
                        << "first is " << first.get_read_number()
                        << ", classes is " << second.get_read_number());
    BOOST_CHECK_MESSAGE(second.get_percentage_useful_reads()
                                == first.get_percentage_useful_reads(),
                        "get_percentage_useful_reads fn don't match."
                        << "first is "
                        << first.get_percentage_useful_reads()
                        << ", classes is "
                        << second.get_percentage_useful_reads());
    BOOST_CHECK_MESSAGE(second.get_fragment_number()
                                == first.get_fragment_number(),
                        "get_fragment_number fn don't match."
                        << "first is " << first.get_fragment_number()
                        << ", classes is "
                        << second.get_fragment_number());
    BOOST_CHECK_MESSAGE(second.get_avg_frag_length()
                                == first.get_avg_frag_length(),
                        "get_avg_frag_length fn don't match."
                        << "first is " << first.get_avg_frag_length()
                        << ", classes is "
                        << second.get_avg_frag_length());
    BOOST_CHECK_MESSAGE(second.get_total_bases()
                                == first.get_total_bases(),
                        "get_total_bases fn don't match."
                        << "first is " << first.get_total_bases()
                        << ", classes is " << second.get_total_bases());
    BOOST_CHECK_MESSAGE(second.get_num_useable_bases()
                                == first.get_num_useable_bases(),
                        "get_num_useable_bases fn don't match."
                        << "first is " << first.get_num_useable_bases()
                        << ", classes is "
                        << second.get_num_useable_bases());
    BOOST_CHECK_MESSAGE(second.get_percentage_useful_bases()
                                == first.get_percentage_useful_bases(),
                        "get_percentage_useful_bases fn don't match."
                        << "first is "
                        << first.get_percentage_useful_bases()
                        << ", classes is "
                        << second.get_percentage_useful_bases());
}

template <class SA> void try_full_alg(std::string file_name,
                  std::unordered_map<std::string, uint_least64_t> input,
                  SA &test_SA) {
    ProblemDescription test_PD_with_file(file_name, "../testdata/Reads.txt", input);

    std::string concat_reads = dump_reads_to_string("../testdata/Reads.txt",
                                            test_PD_with_file.get_SA_length());
    BOOST_TEST(concat_reads.length() == test_SA.get_SA_length());

    std::cout.setstate(std::ios_base::failbit);
    test_SA.run_full_alg();
    std::cout.clear();
    test_order(test_SA, concat_reads);
}

template <class SA> void test_order(const SA &test_SA,
                                       const std::string &concat_reads) {
    if (test_SA.get_SA_length() > 5000000) {
        std::cout << "    Test completed to SA index: " << std::flush;
    }
    for (int i = 0; i < test_SA.get_SA_length()-1; i++) {
        // Assert strings are in the correct order.
        std::string first = get_cyclic_shift_at(test_SA.get_suffix_at(i),
                                                test_SA.get_sort_length(),
                                                concat_reads);
        std::string second = get_cyclic_shift_at(test_SA.get_suffix_at(i+1),
                                                 test_SA.get_sort_length(),
                                                 concat_reads);

        BOOST_REQUIRE_MESSAGE(first.compare(second) <= 0,
                            "i is " << i << " first is " << first
                            <<", " << test_SA.get_suffix_at(i)
                            << " second is " << second << ", "
                            << test_SA.get_suffix_at(i+1));

        if ((test_SA.get_SA_length() > 5000000)
                && !(i%5000000) && (i != 0)) {
            std::cout << i/1000000 << "M  " << std::flush;
        }

        // Assert classes are in the correct order.
        int c_first = test_SA.get_class_at(test_SA.get_suffix_at(i));
        int c_second = test_SA.get_class_at(test_SA.get_suffix_at(i+1));
        BOOST_TEST(c_first <= c_second);
    }
    if (test_SA.get_SA_length() > 5000000) { std::cout << "\n" << std::flush; }
}

int path_inverse(std::stack<bool> path) {
    /*
    * Iterate through a path and calculate the number it produces.
    * Assumed: 1 => add one to the result
    *          0 => double the result.
    */
    int result = 0;

    while (!path.empty()) {
        if (path.top() == 0) {
            result = result*2;
        } else {
            result += 1;
        }
        path.pop();
    }
    return result;
}

std::string get_cyclic_shift_at(int i, int length,
                                const std::string &sequence) {
    std::string result;

    if (i + length - 1 < sequence.length()) {
        // Then the cyclic shift doesn't run over the string length.
        result = sequence.substr(i, length);
    } else {
        result = sequence.substr(i, sequence.length()) +
                 sequence.substr(0, length-(sequence.length()-i));
    }
    return result.substr(0, result.find("."));
}

std::string dump_reads_to_string(std::string file_name, int file_length) {
    std::string test;
    assert(test.max_size() > file_length);
    std::ifstream file(file_name);
    std::string concat_reads((std::istreambuf_iterator<char>(file)),
                              std::istreambuf_iterator<char>() );
    file.close();

    return concat_reads;
}

