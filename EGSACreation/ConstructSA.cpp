// Copyright [2018] Isla Carson
#include <stdint.h>
#include <vector>
#include <bitset>
#include <iterator>
#include <algorithm>
#include <string>
#include <stack>
#include <thread>         // std::thread
#include <mutex>          // std::mutex
#include "ConstructSA.hpp"
std::mutex mtx;
std::mutex mtx_sort_location;           // mutex for critical section

void ConstructSA::validate_user_input(
                std::unordered_map<std::string, uint_least64_t> *input_ptr) {
    /*
    *  Take the user input and validate that it falls within allowable limits.
    *  alg_limit: >= 1  or not given.
    *  num_threads:  1 to 1024 or not given.
    */

    // Give default values for the number of threads and the alg limit.
    input_ptr->emplace("-num_threads", 1);
    input_ptr->emplace("-alg_limit", 80000);

    // Validate alg limit.
    if ((*input_ptr)["-alg_limit"] < 1) {
        std::cout << "-alg_limit must be > 0 or unspecified." << std::endl;
        exit(1);
    }

    // Validate num of threads.
    if ( ((*input_ptr)["-num_threads"] < 1)
         || ((*input_ptr)["-num_threads"] > 1024) ) {
        std::cout << "-num_threads must be > 0 and < 1024 or unspecified."
                  << std::endl;
        exit(1);
    }
}


ConstructSA::ConstructSA(const std::string &fastx,
                         std::unordered_map<std::string, uint_least64_t> input)
            :ProblemDescription(fastx, input, true) {
    /*
    *  Create the ConstructSA class for file fastx and user input.
    *
    *  Load file into Class array.
    *  Get the needed alg steps.
    *  Initialise other variables.
    */

    // Validate that the user input gives a well defined ConstructSA alg.
    std::unordered_map<std::string, uint_least64_t> *input_ptr = &input;
    validate_user_input(input_ptr);

    // Process problem and get the equivalence Classes
    m_classes = process_problem_to_vector();

    if ( get_depth() == 0 ) { set_depth_to_max();
    } else { reduce_depth_to_max(); }
    if (get_precision_req()) { set_precise_alg_steps();
    } else { set_full_alg_steps(); }

    is_valid();
    m_alphabet_size = 5;

    // base indices
    m_base_indices.assign(m_alphabet_size+1, 0);
    m_min_class_for_base.assign(m_alphabet_size, 0);
    find_base_indices();

    // Sort depth and number of equivalence classes initialised.
    m_curr_length = 0;
    num_classes = 0;

    // Suffix Array
    m_SA = std::vector<uint_least64_t> (get_SA_length());
    m_flags = std::vector<bool> (2*m_alphabet_size);

    //  Algorithm Comparison vs CountSort switch point
    m_limit = input["-alg_limit"];

    // Threads
    m_num_threads = input["-num_threads"];
    m_thread_chunks.push_back(0);
    for (int i = 2; i < m_alphabet_size; i++) {
        m_thread_chunks.push_back(m_base_indices[i]);
    }
    m_thread_chunks.push_back(get_SA_length());
    m_threads = std::vector<std::thread> (m_num_threads-1);
    // As the main thread will also do work.
}


// Destructor
ConstructSA::~ConstructSA() {}


void ConstructSA::set_precise_alg_steps() {
    /*
    *  Get the alg steps needed to exactly reach the sort_depth.
    *  0 represents a doubling step, and 1 a +1 step.
    */
    uint_least64_t curr = get_depth();
    std::stack<bool> path;

    while (curr > 1) {
        if (curr%2) {
            // i can't be halved, so +1 needed.
            path.push(1);
            curr = curr - 1;
        } else {
            // Then i is even, so half.
            path.push(0);
            curr = curr/2;
        }
    }

    // End path with '1', as the 1st step is to read and sort for len 1.
    path.push(1);
    m_alg_steps = path;
}


void ConstructSA::set_full_alg_steps() {
    /*
    *  Get the alg steps needed to reach (or exceed) the sort_depth.
    *  0 represents a doubling step, and 1 a +1 step.
    */
    uint_least64_t curr = 1;

    std::stack<bool> path;

    // Double until sort_depth is exceeded.
    while (curr < get_depth()) {
        path.push(0);
        curr *= 2;
    }

    // End path with '1', as the 1st step is to read and sort for len 1.
    path.push(1);
    m_alg_steps = path;
}


std::string ConstructSA::see_alg_steps() const {
    /* Print out alg steps to a string. */
    std::string pretty;

    //  Temp stack to allow iteration.
    std::stack<bool> temp = copy_alg_steps();

    while (!temp.empty()) {
        pretty.append(temp.top()  ? "1" : "0");
        temp.pop();
    }

    return pretty;
}

std::stack<bool> ConstructSA::copy_alg_steps() const {
    /* Get a copy of the alg steps. */
    std::stack<bool> temp = m_alg_steps;
    return temp;
}

/* Public get fns to allow access to class members*/
uint_least64_t ConstructSA::get_sort_length() const {
    return m_curr_length;
}

uint_least64_t ConstructSA::get_class_at(uint_least64_t i) const {
    return m_classes[i];
}

bool ConstructSA::get_flag_at(uint_least64_t i) const {
    return m_flags[i*2];
}


bool ConstructSA::get_sort_choice_for(uint_least64_t i) const {
    return m_flags[i*2+1];
}

bool ConstructSA::get_single_char_flag(uint_least64_t i) const {
    /*
    *  Return the sort-requirement flag for the +1 step.
    *  Only the '.' char with class 0 can be skipped.
    *  All other classes (1 to 4) need sorted.
    *
    *  So: return 1 (no sort) if 0 and 0 otherwise.
    */
    return (i == 0)?1:0;
}

void ConstructSA::set_flag(uint_least64_t i, bool flag) {
    /* Set the sort-requirement flag for a given class. */
    m_flags[i*2] = flag;
}


void ConstructSA::set_sort_choice(uint_least64_t i, bool choice) {
    /*
    * Set the sort-choice flag for a given class.
    * 0 represents comparison sort and 1 represents the combination sort.
    */
    m_flags[i*2+1] = choice;
}

void ConstructSA::find_base_indices() {
    /*
    *  Find the occurrences of each base in the file.
    *  Then find the start of it's suffix interval
    */

    // Iterate through the class vector, counting the occurrence of each char
    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        int index = m_classes[i] + 1;
        if (index < m_alphabet_size) {
                m_base_indices[index] += 1;
        }
    }
    // Get the starting position for each char in the SA for cyclic shifts
    // of length 1.
    for (int i = 1; i < m_alphabet_size; i++) {
        m_base_indices[i] += m_base_indices[i-1];
        m_min_class_for_base[i] = i;
    }
    m_base_indices[m_alphabet_size] = get_SA_length();
}

uint_least64_t ConstructSA::get_base_index(uint_least64_t i) const {
    /* Public function to return the SA start position of each char. */
    return m_base_indices[i];
}

uint_least64_t ConstructSA::get_min_class_for_base(uint_least64_t i) const {
    return m_min_class_for_base[i];
}


uint_least64_t ConstructSA::get_suffix_at(uint_least64_t i) const {
    return m_SA[i];
}

void ConstructSA::update_base_class_start() {
    /*
    *  Find the minimum class for each char.
    *  Used as a comparison to determine the starting char in each class string.
    */
    int num_bases = m_alphabet_size;
    for (int i = 0; i < m_alphabet_size; i++) {
            // Get the class at the start of each base.
            m_min_class_for_base[i] = m_classes[
                                        m_SA[
                                            m_base_indices[i]
                                            ]
                                        ];
    }
}


uint_least64_t ConstructSA::get_immediately_appending_class(uint_least64_t i)
                                                                        const {
    /*
    *   Return the class of the character ( ie the class of the suffix of
    *   length 1 in position i + m_curr_length ) occuring immediately after
    *   the i-th suffix (which has length m_curr_length).
    *
    *   The class can be determined via comparison against the minimum class
    *   for each base.
    *
    *   Note that the return type must be uint_least64_t, in order to match
    *   the return type og get_doubled_appending_class.
    */

    // Written in this manner (as opposed to using %) to control possible
    // outcomes near the MAX of uint_least64_t.
    uint_least64_t appending_suffix = i + m_curr_length;
    appending_suffix = (appending_suffix >= get_SA_length()) ?
                        appending_suffix - get_SA_length() : appending_suffix;
    uint_least64_t appending_class = m_classes[appending_suffix];

    for (int j = 4; j > 0; j--) {
        if (appending_class >= m_min_class_for_base[j]) {
            // Then it is of type i.
            return j;
        }
    }
    // If code gets this far, then class must be 0.
    return 0;
}

uint_least64_t ConstructSA::get_doubled_appending_class(uint_least64_t i) const {
    /*
    *  Get the class of the suffix of length m_curr_length at position i + m_curr_length
    */


    // Written in this manner (as opposed to using %) to control possible
    // outcomes near the MAX of uint_least64_t.
    uint_least64_t appending_suffix = i + m_curr_length;
    appending_suffix = (appending_suffix >= get_SA_length()) ?
                        appending_suffix - get_SA_length(): appending_suffix;
    return m_classes[appending_suffix];
}

int ConstructSA::get_prepending_class(uint_least64_t i) const {
    /*
    *   Get the class of the suffix of length 1 ocurring at position i-1.
    */

    // Written in this manner (as opposed to using %) to control possible
    // outcomes near the MAX of uint_least64_t.
    uint_least64_t prependingsuffix = m_SA[i];
    prependingsuffix = (prependingsuffix == 0) ?
                        prependingsuffix + get_SA_length() - 1 :
                        prependingsuffix - 1;

    uint_least64_t prepending_class = m_classes[prependingsuffix];
    for (int i = 4; i > 0; i--) {
        if (prepending_class >= m_min_class_for_base[i]) {
            // Then it is of type i.
            return i;
        }
    }
    // If code gets this far, then class must be 0.
    return 0;
}



void ConstructSA::initialise_alg() {
    /*
    *  Set up the alg by performing count sort on the char of length 1.
    */

    // Copy indices to a local array so they can be modified without affecting
    // the base count variables.
    std::vector<uint_least64_t> base_index(m_alphabet_size);
    for (int i=0; i < m_alphabet_size; i++) {
        base_index[i] = m_base_indices[i];
        set_flag(i, (i == 0)?1:0);
        set_sort_choice(i, 1);
    }

    // Loop through the reads, sorting the first char of each suffix using a
    // Count Sort. (Note: base incidence already known)
    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        // Get the lexographical standing of the char *i
        int char_standing = m_classes[i];
        // Get it's the suffix location in the SA and assign.
        m_SA[base_index[char_standing]] = i;
        // Increment
        base_index[char_standing] +=1;
    }

    // Now sorted to length 1.
    m_curr_length = 1;
    num_classes = m_alphabet_size;
}


void ConstructSA::update_SA_comb_sort(uint_least64_t class_start,
                uint_least64_t class_end,
                std::unordered_map<uint_least64_t, uint_least64_t> &count_map,
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const) {
    /*
    *   Use the incidence of classes (stored in count_map) to comb sort the
    *   section of the SA between class_start and class_end.
    */

    if (count_map.size() > 1) {
        // Order the classes and cummulatively total to get the start index of
        // their suffix intervals.
        order_and_total(count_map);

        // Create a placeholder vector of the required length.
        std::vector<uint_least64_t> new_SA(class_end - class_start);

        // Perform a count sort.
        for (uint_least64_t i = class_start; i < class_end; i++) {
            uint_least64_t appending_suffix_class = (this->*get_class)(m_SA[i]);
            uint_least64_t new_index = count_map[appending_suffix_class];
            new_SA[new_index] = m_SA[i];

            count_map[appending_suffix_class] += 1;
        }

        // Copy the completed new_SA into m_SA.
        for (uint_least64_t i = class_start; i < class_end; i++) {
            m_SA[i] = new_SA[i-class_start];
        }
    }
}

void ConstructSA::order_and_total(
                std::unordered_map<uint_least64_t, uint_least64_t> &count_map) {
    /*
    *   Order the classes and cummulatively total to get the start index of
    *   their suffix intervals.
    */

    // Extract the keys from the count map
    std::vector<uint_least64_t> ordered_keys(count_map.size());
    uint_least64_t i = 0;
    for (auto it = count_map.begin(); it != count_map.end(); ++it) {
        ordered_keys[i] = it->first;
        i++;
    }

    // Order the keys via a comparison sort. O(num_keys * log(num_keys))
    std::sort(ordered_keys.begin(), ordered_keys.end());

    // Iterate over map to cumulatively total over the ordered keys.
    uint_least64_t running_total = 0;
    for (uint_least64_t i = 0; i < count_map.size(); i++) {
        uint_least64_t holder = count_map[ordered_keys[i]];
        count_map[ordered_keys[i]] = running_total;
        running_total += holder;
    }
}

int ConstructSA::get_next_work_chunk(std::vector<bool> *sort_commenced) {
    /*
    *  Get the index of the next chunk of work, and mark it as started.
    *  Return -1 if no work chunks remaining.
    */
    mtx_sort_location.lock();
    for (int i = 0; i < sort_commenced->size(); i++) {
        if (!(sort_commenced->at(i))) {
            sort_commenced->at(i) = true;
            mtx_sort_location.unlock();
            return i;
        }
    }
    mtx_sort_location.unlock();
    return -1;
}

uint_least64_t ConstructSA::get_class_end(uint_least64_t i) {
    /*  Find the end of the current class within the SA. */

    uint_least64_t current_class = m_classes[m_SA[i]];
    while ((i < get_SA_length()) && (m_classes[m_SA[i]] == current_class)) {
        // Increase class_end until the current class is done.
        i++;
    }
    return i;
}

uint_least64_t ConstructSA::count_and_get_class_end(uint_least64_t i,
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                std::unordered_map<uint_least64_t, uint_least64_t> &count_map) {
    /*
    *   Find the end of the current class within the SA, whilst counting the
    *   occurences of each appending class.
    */

    uint_least64_t current_class = m_classes[m_SA[i]];
    while ((i < get_SA_length()) &&  (m_classes[m_SA[i]] == current_class)) {
        // Add 1 to the appending_suffix count.
        uint_least64_t appending_suffix = (this->*get_class)(m_SA[i]);
        count_map.emplace(appending_suffix, 0);
        count_map[appending_suffix] += 1;

        // Increase class_end until the current class is done.
        i++;
    }
    return i;
}

void ConstructSA::sort_suffix_pairs(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const) {
    /*
    *   Sort the SA of suffices of length 2*m_curr_length or 1+m_curr_length
    *   (alg_step type dependent).  Allocate each chunk of work to a thread.
    */

    // A vector to record which chunks of works have been started.
    std::vector<bool> sort_commenced;
    for (int i = 0; i < m_thread_chunks.size()-1; i++) {
        sort_commenced.push_back(false);
    }

    // Number of threads needed.
    int threads = (m_num_threads < m_thread_chunks.size()-1) ?
                   m_num_threads : m_thread_chunks.size()-1;

    // Spawn threads
    for (int j = 0; j < threads-1; ++j) {
        m_threads[j] = std::thread(&ConstructSA::sort_thread, this,
                                   get_class, &sort_commenced);
    }
    // Give the main thread work to do too.
    sort_thread(get_class, &sort_commenced);

    // Join all threads & don't preceed until they're all done.
    for (int j = 0; j < threads-1; ++j) { m_threads[j].join(); }
}


void ConstructSA::sort_thread(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                std::vector<bool> *sort_pointer) {
    /*
    *   Manage the thread as it works on chunks of m_SA.
    *   The thread should:
    *       -- Fetch a chunk of work (ie an interval of m_SA).
    *       -- Iterate through the interval, identifying the class boundaries.
    *       -- For each sub-interval (as defined by the class boundaries),
    *          dispatch to the relevant sort function.
    *       -- Sort fns are called within the thread.
    */

    // Get the next available chunk of work.
    int work_chunk = get_next_work_chunk(sort_pointer);

    while (work_chunk != -1) {
        // Get the boundaries of the work chunk in m_SA.
        uint_least64_t class_start = m_thread_chunks[work_chunk];
        uint_least64_t work_chunk_end = m_thread_chunks[work_chunk + 1];

        // Within the work chunk, there's multiple subintervals.
        // Cycle through them until the end of the work chunk is reached.
        while (class_start < work_chunk_end) {
            uint_least64_t current_class = m_classes[m_SA[class_start]];
            uint_least64_t class_end;

            if (!get_flag_at(current_class)) {
                // Sorting is required.

                // Find the sort alg choice for the current class.
                bool alg_choice = get_sort_choice_for(current_class);

                if (alg_choice) {
                    // Get class_end and counts of the appending suffices.
                    std::unordered_map<uint_least64_t, uint_least64_t> count_map;
                    class_end = count_and_get_class_end(class_start, get_class,
                                                        count_map);
                    // Perform comb sort.
                    update_SA_comb_sort(class_start, class_end,
                                        count_map, get_class);
                } else {
                    // Get the end of the class
                    class_end = get_class_end(class_start);

                    // Get iterator pointers to the m_SA interval to be sorted.
                    std::vector<uint_least64_t>::iterator start = m_SA.begin();
                    std::advance(start, class_start);
                    std::vector<uint_least64_t>::iterator end = m_SA.begin();
                    std::advance(end, class_end);

                    // Sort using comparison sort.
                    std::sort(start, end,
                        [this, get_class](uint_least64_t i, uint_least64_t k) {
                            return (this->*get_class)(i)
                                                    < (this->*get_class)(k);
                        });
                }
            } else {
                // No sort required.  So skip to the end of the class.
                class_end = get_class_end(class_start);
            }
            class_start = class_end;
        }

        // Get the next available chunk of work.
        work_chunk = get_next_work_chunk(sort_pointer);
    }
}

void ConstructSA::revise_classes_thread(
                        std::vector<std::vector<uint_least64_t>> *class_section,
                        std::vector<uint_least64_t> *counter,
                        std::vector<bool> *work_pointer) {
    /*
    *  Thread to update m_classes.  Thread should:
    *       Get the next work chunk.
    *       Iterate through the section of m_classes between
            class_section[work_chunk][0] and class_section[work_chunk][len-1],
            updating each entry to contain the needed class.
    */

    // Get the next available chunk of work.
    int work_chunk = get_next_work_chunk(work_pointer);

    while (work_chunk != -1) {
        // Get the boundaries of the work chunk in m_SA.
        for (uint_least64_t j = 0;
             j < class_section->at(work_chunk).size()-1;
             j++) {
            // Iterate through each interval in the class section.
            for (uint_least64_t k = class_section->at(work_chunk)[j];
                 k < class_section->at(work_chunk)[j+1];
                 k++) {
                // Update each class to j +num  of previous classes.
                m_classes[m_SA[k]] = counter->at(work_chunk) + j;
            }
        }
        // Get the next work chunk.
        work_chunk = get_next_work_chunk(work_pointer);
    }
}


void ConstructSA::update_classes(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                bool (ConstructSA::*get_flag)(uint_least64_t) const,
                bool mode) {
    /*
    *   Manage several threading processes in order to update the classes of
    *   the suffices.
    */

    // Set up threads.
        // Vector to record which chunks have been worked on.
    std::vector<bool> work_commenced;
    int num_work_chunks = m_thread_chunks.size()-1;
    for (int i = 0; i < num_work_chunks; i++) {
            work_commenced.push_back(false);
        }
        // Vectors to record the results of each thread.
    std::vector<std::vector<uint_least64_t>> class_section(num_work_chunks);
    std::vector<std::vector<bool>> flag_section(num_work_chunks);

    // Get the number of threads.
    int threads = (m_num_threads < num_work_chunks) ?
                                            m_num_threads : num_work_chunks;
    // Spawn the threads.
    for (int j = 0; j < threads-1; j++) {
        // Get class changes.
        m_threads[j] = std::thread(&ConstructSA::get_class_thread, this,
                                   get_class, get_flag,
                                   &work_commenced, &class_section,
                                   &flag_section);
    }
    // Give the main thread work to do too.
    get_class_thread(get_class, get_flag,
                     &work_commenced, &class_section, &flag_section);
    // Join the threads.
    for (int j = 0; j < threads-1; ++j) { m_threads[j].join(); }

    // Merge flag_sections into one vector, m_flags - the class member.
    m_flags.clear();
    uint_least64_t total_flags = 0;
    for (int i = 0; i < num_work_chunks; i++) {
        total_flags += flag_section.size();
    }
    m_flags.reserve(total_flags);
    for (int i = 0; i < num_work_chunks; i++) {
        std::move(flag_section[i].begin(),
                  flag_section[i].end(),
                  std::back_inserter(m_flags));
    }

    // Process class_sections into m_classes.
    std::vector<uint_least64_t> counter(num_work_chunks);
    counter[0] = 0;
    for (int i = 1; i < num_work_chunks; i++) {
        counter[i] = counter[i-1] + class_section[i-1].size()-1;
    }

    // Reset the work pointer.
    for (int i = 0; i < num_work_chunks; i++) { work_commenced[i] = false; }
    // Spawn threads
    for (int j = 0; j < threads-1; j++) {
        // Revise m_classes to record the new classes.
        m_threads[j] = std::thread(&ConstructSA::revise_classes_thread, this,
                                   &class_section, &counter, &work_commenced);
        }
    // Give the main thread work to do too.
    revise_classes_thread(&class_section, &counter, &work_commenced);
    // Join the threads.
    for (int j = 0; j < threads-1; ++j) { m_threads[j].join(); }

    // Update the record of the minimum class value for each base.
    update_base_class_start();

    // Print stats out to the user.
    double avg_classes = num_classes;
    num_classes = m_classes[m_SA[get_SA_length()-1]];
    std::cout << "Number of classes: " << num_classes
              << "\nAvg number of unique app_class per prev class:"
              << static_cast<double>(num_classes)/avg_classes
              << "\nAvg number of suffices per current class: "
              << static_cast<double>(get_SA_length())/
                    static_cast<double>(num_classes)
              << "\nAvg new class size: "
              << static_cast<double>(get_SA_length())/
                    static_cast<double>(num_classes)
              << "\n" << std::flush;


    // Update the current length.
    if (mode) {
        // Then step is append by 1. So:
        m_curr_length +=1;
    } else {
        // Then step is a doubling of suffix length.
        m_curr_length *= 2;
    }

    // Use the class_sections vector to reduce the size of the work chunks.
    if (m_curr_length < 8) {
        // If a class directly follows one that requires no further sorting,
        // it doesn't require an entry. This is to eliminate work chunks
        // consisting of only classes that don't required further sorting.
        bool requires_work_chunk = 1;

        m_thread_chunks.clear();
        for (int i = 0; i < num_work_chunks; i++) {
            for (uint_least64_t j = 0; j < class_section[i].size()-1; j++) {
                if (requires_work_chunk) {
                    m_thread_chunks.push_back(class_section[i][j]); }

                // Update requires_work_chunk to indicate if the next class
                // needs it a new work chunk.
                if (get_flag_at(counter[i]+j)) { requires_work_chunk = 0;
                } else { requires_work_chunk = 1; }
            }
        }
        m_thread_chunks.push_back(get_SA_length());
    }
}


uint_least64_t ConstructSA::get_classes_in_interval(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                bool (ConstructSA::*get_flag)(uint_least64_t) const,
                uint_least64_t class_start,
                std::vector<uint_least64_t> &class_changes,
                std::vector<bool> &flag_changes) {
    /*
    *   Get the class changes within an interval where the first of the pair
    *   of suffices is set.  When the first suffix in the pair is found to
    *   have a different class from the previous first suffix in the pair, the
    *   function ends and returns the index of the suffix with the new class.
    *
    *   AKA return the end of the current class.
    */

    // Info on current suffix.
    uint_least64_t current_class = m_classes[m_SA[class_start]];
    bool current_flag = get_flag_at(current_class);
    uint_least64_t app_class = (this->*get_class)(m_SA[class_start]);
    bool app_flag = (this->*get_flag)(app_class);  // dependent on mode

    // Add info for the current suffix to the class and flag changes.
    class_changes.push_back(class_start);
    flag_changes.push_back(app_flag || current_flag);

    // The class currently contains 1 member.
    uint_least64_t counter = 1;

    // Setup for inner loop. (class_start suffix already dealt with).
    uint_least64_t i = class_start + 1;
    uint_least64_t prev_app_class = app_class;  // now to be compared against.

    // Find the end index of the current class.
    while ((i < get_SA_length()) &&
          (m_classes[m_SA[i]] == current_class)) {
        // In this loop, the class of the start of the suffix is the same.
        // So only the second suffixes need to be compared.
        uint_least64_t curr_app_class = (this->*get_class)(m_SA[i]);

        // Compare this appending suffix class against the prev.
        if (!current_flag  &&
            ((prev_app_class != curr_app_class) || class_start == i)) {
            // Got the length of the previously occurring class.
            flag_changes.push_back((counter <= m_limit)?0:1);
            // Reset the counter
            counter = 0;

            // New Class
            class_changes.push_back(i);
            app_flag = (this->*get_flag)(curr_app_class);
            // Push back flags.
            flag_changes.push_back(app_flag || current_flag);

            // update the prev class for comparison.
            prev_app_class = curr_app_class;
        }

        // Increase the class counter.
        counter += 1;

        // Increase i until the current class is done.
        i++;
    }

    // Got the length of the previously occurring class.
    flag_changes.push_back((counter <= m_limit)?0:1);

    return i;
}

void ConstructSA::get_class_thread(
                uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const,
                bool (ConstructSA::*get_flag)(uint_least64_t) const,
                std::vector<bool> *work_pointer,
                std::vector<std::vector<uint_least64_t>> *class_section_ptr,
                std::vector<std::vector<bool>> *flag_section_ptr) {
    /*
    *   Manage a thread as it finds the class changes within each work chunk.
    *   The thread should:
    *       Get the next work chunk
    *       Iterate through the work chunk, finding the start of each current
    *       class.
    *       If no further sorting is required, record the class and continue.
    *       Otherwise, send the class to the get_classes_in_interval fn (which
    *       returns the end of
    *       current class).
    *       Continue until the end of the work chunk is reached, and then
    *       fetch the next chunk.
    */

    // Get the next available chunk of work.
    int work_chunk = get_next_work_chunk(work_pointer);

    while (work_chunk != -1) {
        // Get the boundaries of the work chunk in m_SA.
        uint_least64_t class_start = m_thread_chunks[work_chunk];
        uint_least64_t work_chunk_end = m_thread_chunks[work_chunk + 1];

        // Within the work chunk, there's multiple subintervals.
        // Cycle through them until the end of the work chunk is reached.
        while (class_start < work_chunk_end) {
            uint_least64_t current_class = m_classes[m_SA[class_start]];
            uint_least64_t class_end;

            if (!get_flag_at(current_class)) {
                // Further sorting will be required.
                // Get the class end and update the list class change locations.
                class_end = get_classes_in_interval(
                                            get_class, get_flag, class_start,
                                            class_section_ptr->at(work_chunk),
                                            flag_section_ptr->at(work_chunk));
            } else {
                // No further sorting required, so do nothing but record that
                // this class exists:
                class_section_ptr->at(work_chunk).push_back(class_start);
                // Record it's flags.
                flag_section_ptr->at(work_chunk).push_back(1);
                flag_section_ptr->at(work_chunk).push_back(0);
                // And then skip to the end of the class.
                class_end = get_class_end(class_start);
            }
            class_start = class_end;
        }
        class_section_ptr->at(work_chunk).push_back(work_chunk_end);

        // Get the next work chunk.
        work_chunk = get_next_work_chunk(work_pointer);
    }
}


void ConstructSA::run_full_alg() {
    /*
    *   Run the full algorithm, by iterating through alg_steps, selecting the
    *   class and flag functions for the given mode, calling sort_suffix_pairs
    *   and then update_classes.
    */

    // Get the alg steps
    std::stack<bool> alg_path = copy_alg_steps();

    // Initialise alg, and first sort from the path.
    initialise_alg();
    alg_path.pop();
    std::cout << "\nSearch Completed to 1\n";

    while (!alg_path.empty()) {

        // Loop through the path stack, performing the neccessary sort mode.
        uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const;
        bool (ConstructSA::*get_flag)(uint_least64_t) const;
        bool mode = alg_path.top();

        // Define get_class function dependant on mode.
        if (mode) {
            // true, so +1
            get_class = &ConstructSA::get_immediately_appending_class;
            get_flag = &ConstructSA::get_single_char_flag;
        } else {
            // false, so x2
            get_class = &ConstructSA::get_doubled_appending_class;
            get_flag = &ConstructSA::get_flag_at;
        }

        sort_suffix_pairs(get_class);

        update_classes(get_class, get_flag, mode);
        alg_path.pop();
        std::cout << "Completed to " <<m_curr_length << "\n--------\n"
                  << std::flush;
    }
}



void ConstructSA::sort_suffix_pairs(bool mode) {
    /* Interface to sort_suffix_pairs for testing when only mode is known. */
    uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const;
    if (mode) {
            // true, so +1
            get_class = &ConstructSA::get_immediately_appending_class;
        } else {
            // false, so x2
            get_class = &ConstructSA::get_doubled_appending_class;
        }
    sort_suffix_pairs(get_class);
}

void ConstructSA::update_classes(bool mode) {
    /* Interface to update_classes for testing when only mode is known. */
    uint_least64_t (ConstructSA::*get_class)(uint_least64_t) const;
    bool (ConstructSA::*get_flag)(uint_least64_t) const;
    if (mode) {
            // true, so +1
            get_class = &ConstructSA::get_immediately_appending_class;
            get_flag = &ConstructSA::get_single_char_flag;
        } else {
            // false, so x2
            get_class = &ConstructSA::get_doubled_appending_class;
            get_flag = &ConstructSA::get_flag_at;
        }
    update_classes(get_class, get_flag, mode);
}


void ConstructSA::set_all_sort_choices(bool choice) {
    /*
    *   Change all sort alg choices to be 'choice'.
    *   This allows for testing of each alg separately.
    */

    for (uint_least64_t i =0; i < (m_flags.size())/2; i++) {
        set_sort_choice(i, choice);
    }
}

void ConstructSA::create_flagged_fm_index(std::vector<bool> *lcp_flags,
                                    std::vector<uint_least64_t> *fm_index) {
    /*
    *  Convert the current
    */

    lcp_flags->reserve(get_SA_length());

    // Assert precise-depth sorting is used.
    if (!get_precision_req()) {
        std::cerr << "The lcp flagged bwt cannot be created when imprecise "
                  << "sort depth is used.  Please sort precisely or create "
                  << "only the BWT.";
        exit(0);
    }

    // Each class shares a prefix up to a length d or the first $.
    // If the class shares a prefix smaller than d (i.e. it matches only up to the first $), the alg won't use that SA for patttern length reduction.  Hence, it's 'safe' to count it as having a full match to d.

    //  For each class, set the LCP array to be either 0 or 1 for the full length.  Then sort the next char, so that patterns of length d+1 can be searched for.

    // First, run through m_SA, alternating LCP flag at the end of each class.
    bool flag = 0;
    int counter = 0;
    int avg_counter = 0;
    int num_classes_counter = 0;
    uint_least64_t prev_class = 0;
    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        uint_least64_t current_class = m_classes[m_SA[i]];
        if (current_class != prev_class) {
            flag = !flag;
            // if ( get_flag_at(current_class) && !get_flag_at(prev_class) ) {
            //     // Then the termination char has been passed.
            //     num_classes_counter++;
            // }
        }
        // if ( get_flag_at(current_class) ) {
        //     // Then the termination char has been passed.
        //     counter++;
        //     avg_counter++;
        // }

        (*lcp_flags)[i] = flag;
        prev_class = current_class;
    }


    // std::cout << "\n\n" << "$ PASSED COUNTER IS AT " << counter << " OF "
    //           << get_SA_length()
    //           << "\n There are " << num_classes_counter
    //           << " such classes of "
    //           << m_classes[m_SA[get_SA_length()-1]] << " total. "
    //           << "\n The average size of these classes is "
    //           << static_cast<double>(avg_counter)
    //                         /static_cast<double>(num_classes_counter)
    //           << "\n\n";

    std::cout << "LCP flags written." << "\n--------\n" <<  std::flush;

    // Second, perform a sort of depth + 1 on the remaining classes.
    sort_suffix_pairs(true);
    update_classes(true);
    std::cout << "Sort completed to " << m_curr_length
              << std::flush;

    // Create the FM-index from the new-ly sorted classes.
    create_fm_index(fm_index);

    std::cout << "\n--------\n" << "FM Index created."
              << "\n--------\n" <<  std::flush;
}


void ConstructSA::create_fm_index(std::vector<uint_least64_t> *fm_index) {
    /*
    *  Convert the current
    */

    fm_index->reserve(get_SA_length()*5);

    // All counts are initially 0.
    for (int i = 0; i < 5; i++) { (*fm_index)[i] = 0; }

    for (uint_least64_t j = 0; j < get_SA_length(); j++) {
        // Initialise each row with the previous.
        for (int i = 0; i < 5; i++) {
            (*fm_index)[5*(j+1)+i] = (*fm_index)[5*j+i];}
        // Add the current char to the count.
        (*fm_index)[5*(j+1)+get_prepending_class(j)] += 1;
    }
}

std::string ConstructSA::create_BWT() {
    std::string bwt(get_SA_length(), '.');
    for (uint_least64_t j = 0; j < get_SA_length(); j++) {
        bwt[j] = convert_class_to_base(get_prepending_class(j));
    }
    return bwt;
}

void ConstructSA::print_flagged_fm_index(const std::string SA_file) {
    std::vector<bool> lcp_flags;
    std::vector<uint_least64_t> fm_index;

    std::ofstream write_output(SA_file);
    if (!write_output.is_open()) {
        std::cerr << "SA output file cannot be opened." << std::endl;
        exit(1);
    }
    write_output << "Meta-Info\n";
    write_output << get_fastx_file_name() << "---read origin file\n";
    write_output << m_curr_length << "---sort depth of lcp flags\n";
    write_output << get_SA_length() << "---suffix array length\n";
    write_output << get_fragment_number() << "---number of fragments\n";
    write_output << num_classes
                 << "---number of classes at lcp flag sort depth\n";

    write_output << "-------------\nlcp_flag,SA";
    create_flagged_fm_index(&lcp_flags, &fm_index);
    // std::string bwt = create_BWT();


    // write_output << "," << bwt[0];
    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        write_output << "\n" << lcp_flags[i];
        for (int j = 0; j < 5; j++) {
            write_output << "," << fm_index[5*i + j];
        }
        // write_output << "," << bwt[i];
    }

    // Last row - No LCP nor BWT
    write_output << "\n" << !lcp_flags[get_SA_length()-1];
    for (int j = 0; j < 5; j++) {
        write_output << "," << fm_index[5*(get_SA_length()) + j];
    }


    write_output.close();
}



void ConstructSA::print_fm_index(const std::string SA_file) {
    std::vector<uint_least64_t> fm_index;

    std::ofstream write_output(SA_file);
    if (!write_output.is_open()) {
        std::cerr << "SA output file cannot be opened." << std::endl;
        exit(1);
    }

    write_output << "SA";
    create_fm_index(&fm_index);

    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        write_output << "\n";
        for (int j = 0; j < 5; j++) {
            write_output << "," << fm_index[5*i + j];
        }
    }

    write_output.close();
}


uint_least64_t ConstructSA::get_compressed_size() {

    // At a length = kmer size, record if term char passed
    m_term_char_passed.reserve(get_SA_length());
    for (uint_least64_t i = 0; i < get_SA_length(); i++) {
        int current_class = m_classes[i];
        if (get_flag_at(current_class)) {
            // The term char has been passed
            m_term_char_passed[i] = 1;
        } else {
            m_term_char_passed[i] = 0;
        }
    }
    // Sort once more, for k+1
    sort_suffix_pairs(true);
    update_classes(true);
    std::cout << "Sort completed to " << m_curr_length
              << std::flush;


    // Get compression stats
    uint_least64_t compressed_length_no_wee_kmers = 0;
    uint_least64_t compressed_length = 0;
    uint_least64_t old_compressed_length = 0;
    uint_least64_t i = 0;

    int current_base = -1;
    int prev_base = -1;
    while (i < get_SA_length()) {
    // RLE
        current_base = get_prepending_class(i);
        if (current_base != prev_base) {old_compressed_length++;
        }
        prev_base = current_base;
        i++;
    }

    std::vector<bool> char_no_wee_kmers = {0, 0, 0, 0, 0};

    int current_class;
    i=0;
    while ( i < get_SA_length() ) {
        // Get next class
        current_class = m_classes[m_SA[i]];

        // Reset bools.
        std::vector<bool> char_representation = {0, 0, 0, 0, 0};

        // Iterate through class.
        int current_base;
        while (m_classes[m_SA[i]] == current_class) {
            current_base = get_prepending_class(i);
            char_representation[current_base] = 1;
            char_no_wee_kmers[current_base] = 1;
            i++;
        }

        bool curr_flag = m_term_char_passed[m_SA[i-1]];
        bool next_flag=0;
        if (i<get_SA_length()) {
            next_flag = m_term_char_passed[m_SA[i]];
        }
        // Add the number of represented chars in the BWT to the length.
        for (int i = 0; i < 5; i++) {
            compressed_length += char_representation[i];
            if ( !curr_flag || (curr_flag && !next_flag)) {
                // if !curr_flag, then the term flag hasn't been passed, so record.  If curr_flag, then term char not passed.  So record iff the next_flag differs.  I.e. if curr_flag=1 and next_flag=0
                compressed_length_no_wee_kmers
                                += char_representation[i];
                // clear bools
                char_no_wee_kmers[i] = 0;
            }
        }
    }

    std::cout << "\n Compressed BWT length is "
              << old_compressed_length << " of " << get_SA_length()
              << " total, or "
              << (static_cast<double>(old_compressed_length)/
                        static_cast<double>(get_SA_length()))*100
              << "%.\n";

    std::cout << "\n Class Bucket Compressed BWT length is "
              << compressed_length << " of " << get_SA_length()
              << " total, or "
              << (static_cast<double>(compressed_length)/
                        static_cast<double>(get_SA_length()))*100
              << "%.\n";

    std::cout << "\n No wee k-mers, class Bucket Compressed BWT length is "
              << compressed_length_no_wee_kmers << " of " << get_SA_length()
              << " total, or "
              << (static_cast<double>(compressed_length_no_wee_kmers)/
                        static_cast<double>(get_SA_length()))*100
              << "%.\n";

    return compressed_length;
}




int ConstructSA::convert_base_to_class(char base) const {
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

char ConstructSA::convert_class_to_base(int char_eq) const {
    switch (char_eq) {
        case 1: return 'a';
        case 2: return 'c';
        case 3: return 'g';
        case 4: return 't';
        case 0: return '.';
        default: return 'x';
    }
}

