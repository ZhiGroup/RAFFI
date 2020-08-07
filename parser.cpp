/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for parsing the output of RaPID.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <boost/thread.hpp>
#include <future>
#include <mutex>
#include <unordered_map>
#include <condition_variable>

#include <stdio.h>
#include <string.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "mapper.hpp"
#include "parser.hpp"
#include "classifier.hpp"
#include "proceed.hpp"
#include "dumpable.hpp"
#include "ordering.hpp"
#include "RaPIDaffin.hpp"
#include <vector>

// Number of individuals to process on one chromosome between two synchronizations
#define NUM_IDS_PER_CYCLE 1000

// Information parsed from a line.
struct line_info {
    int id1_index = -1;
    int id2_index = -1;
    int hap1;
    int hap2;
    int starting_site;
    int ending_site;
};

static inline int haps_to_encoding(int hap1, int hap2);
static inline int haps_encoding_to_complement(int encoding);
static inline bool intersect(int intersection_start, int intersection_end);
static inline int get_intersection_start(int start1, int start2);
static inline int get_intersection_end(int end1, int end2);
static inline std::unique_ptr<std::vector<std::pair<int, int>>> merge_two_segment_vectors(
    std::vector<std::pair<int, int>> &segments1,
    std::vector<std::pair<int, int>> &segments2);
static inline std::unique_ptr<std::vector<std::pair<int, int>>> merge_four_segment_vectors(
    std::vector<std::pair<int, int>> &segments1,
    std::vector<std::pair<int, int>> &segments2,
    std::vector<std::pair<int, int>> &segments3,
    std::vector<std::pair<int, int>> &segments4);
static inline double compute_total_ibd1(
    std::vector<std::pair<int, int>> &segments,
    int chromosome_number);
static void parse_line(
    std::string &line,
    struct line_info &info,
    class Ordering &order);
static void update_total_ibd1(
    int chromosome_number,
    int id_index,
    std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<int, int>>>> &id_to_haps_to_segment,
    std::unordered_map<int, std::unordered_map<int, struct pair_stats>> &matrix);
static bool process_segment(
    struct line_info &info,
    int chromosome_number,
    int chromosome_start,
    std::vector<int> &prev_ids,
    std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<int, int>>>> &id_to_haps_to_segment,
    std::unordered_map<int, std::unordered_map<int, struct pair_stats>> &matrix);
static void worker(
    int thread,
    int num_threads,
    Proceed &proceed,
    Dumpable &dumpable,
    class Ordering &order,
    std::string &rapid_output_path,
    int chromosome_start,
    int chromosome_end,
    std::unordered_map<int, std::unordered_map<int, struct pair_stats>> &matrix);
static int get_min_kinship_coefficient(int max_degree);

/**
 * Master thread responsible for synchronization, writing to either temporary or final
 * output at synchronization, and reading back and inferring pairs written to temporary
 * output.
 *
 * @param vcf_path path to a VCF file.
 * @param rapid_output_path folder that stores the outputs of RaPID.
 *     Assume output of Chromosome i is stored in subfolder i.
 *     E.g. output of Chromosome 1 stored in {rapid_output_path}/1/
 * @param map_path folder that stores genetic maps.
 *     Assume the map for Chromosome i is named as chr{i}.rMap.
 *     E.g. chr22.rMap for chromosome 22.
 * @param max_degree largest degree user is looking for
 * @param num_threads number of worker threads to spawn.
 *     Each thread is responsible for 22 / num_threads chromosomes.
 * @param out final output
 */
void
master(
    std::string &vcf_path,
    std::string &rapid_output_path,
    std::string &map_path,
    int max_degree,
    unsigned int num_threads,
    std::ostream &out)
{
    output_header(out);

    // Initialize genetic maps
    init_maps(map_path);

    // Maximum number of threads is 22
    num_threads = std::min(num_threads, (unsigned int) NUM_CHROMOSOMES);
    int num_chromosomes_per_thread = NUM_CHROMOSOMES / num_threads;

    std::vector<std::future<void>> futures;
    futures.reserve(num_threads);

    Proceed proceed(num_threads);
    Ordering id_ordering(vcf_path);
    Dumpable dumpable_index(id_ordering);

    // An vector of matrix. One for each thread. Each matrix is a
    // 2D unordered map where M[i][j] is a struct pair_stats that
    // records the total IBD1 and IBD2 between individual
    // with index i and individual with index j.
    std::vector<std::unordered_map<int, std::unordered_map<int, struct pair_stats>>> matrices(num_threads);

    uint64_t num_dumped = 0;
    {
        // Temporary output
        std::ofstream temp_file_out = std::ofstream(".temporary", std::ios::out | std::ios::trunc | std::ios::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::output> buffer_out;
        buffer_out.push(boost::iostreams::gzip_compressor());
        buffer_out.push(temp_file_out);
        std::ostream temp_out(&buffer_out);

        // Pairs with kinship coefficents smaller than this will not be written to any output.
        double min_kinship_coefficient = get_min_kinship_coefficient(max_degree);

        // Spwan worker threads
        for (unsigned int thread = 0; thread < num_threads; ++thread) {
            futures.push_back(std::async(
                std::launch::async,
                worker,
                thread,
                num_threads,
                std::ref(proceed),
                std::ref(dumpable_index),
                std::ref(id_ordering),
                std::ref(rapid_output_path),
                thread * num_chromosomes_per_thread + 1,
                // Last thread handles all remaining chromosomes
                thread == num_threads - 1 ? NUM_CHROMOSOMES : (thread + 1) * num_chromosomes_per_thread,
                std::ref(matrices[thread])
            ));
        }

        int count = 0;
        bool done = false;
        while (!done) {
            std::unique_lock<std::mutex> lock(proceed.mutex);

            // Wait until all threads are blocked
            while (!proceed.has_all_blocked()) {
                proceed.wait_for_workers(lock);
            }

            // Write dumpable individuals to output
            std::pair<int, int> range = dumpable_index.get_dumpable_indices();
            num_dumped += dump_range(
                max_degree,
                min_kinship_coefficient,
                id_ordering,
                range,
                matrices,
                temp_out,
                out
            );
            count += range.second - range.first + 1;

            std::cout << count << " individuals processed\r";
            std::cout.flush();

            // Update index of last dumpable individual
            dumpable_index.set_previous_last_dumpable_index(range.second);

            // Reset number of blocked threads to number of unfinished threads
            proceed.update_num_blocked();
            proceed.allow_all_threads_proceed();

            done = proceed.has_all_finished();

            // Resume worker threads
            proceed.signal_workers();
        }
    }

    // Free genetic maps
    deinit_maps();

    std::cout << std::endl;
    std::cout << "Wrote " << num_dumped << " candidate pairs to disk";

    // Wait for all threads to finish
    for (std::future<void> &f : futures) {
        f.get();
    }

    {
        // Open temporary file for reading
        std::ifstream temp_file_in(".temporary", std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> buffer_in;
        buffer_in.push(boost::iostreams::gzip_decompressor());
        buffer_in.push(temp_file_in);
        std::istream temp_in(&buffer_in);

        // Adjust inference boundaries
        shift_boundary();

        // Read in candidate pairs and infer relatedness based on adjusted boundaries
        infer_candidates(max_degree, num_dumped, temp_in, id_ordering, out);
    }

    // Remove temporary file
    if (std::remove(".temporary")) {
        throw std::runtime_error {"Failed to remove temporary file"};
    }

}


/**
 * @param max_degree largest degree user is looking for
 *
 * @return minimum kinship coefficient a pair should reach to be at least max_degree
 *     based on MIN_POWER, the minimum acceptable power of RaPID.
 */
static int get_min_kinship_coefficient(int max_degree) {
    double min_kinship_coefficient;
    switch (max_degree) {
        case 1:
            min_kinship_coefficient = PO_FS_START * MIN_POWER;
            break;
        case 2:
            min_kinship_coefficient = SECOND_START * MIN_POWER;
            break;
        case 3:
            min_kinship_coefficient = THIRD_START * MIN_POWER;
            break;
        case 4:
            min_kinship_coefficient = FOURTH_START * MIN_POWER;
            break;
        default:
            throw std::runtime_error {"Degrees less than 1 or beyond 4 are not supported"};
    }
    return min_kinship_coefficient;
}


/**
 * Worker thread responsible for parsing one or more chromosomes.

 * For each chromosome, parse NUM_IDS_PER_CYCLE individuals between two synchronizations.
 * Block after this is done for all the chromosomes this thread is responsible for.
 * Last thread that reaches synchronization (blocked) should resume the master thread.
 * Restart after master thread writes dumpable individuals to output.
 * Repeat the cycle until the chromosomes have been exausted.
 * Terminate after the chromosomes have been exausted.
 *
 * @param thread thread ID.
 * @param num_threads total number of worker threads.
 * @param proceed a Proceed that determines if this worker thread can resume.
 * @param dumpable a Dumpable that determines ranges of indices of individuals
 *     that can be written to output.
 * @param order an Ordering that specifies the ordering of the IDs as they appear in VCF.
 * @param rapid_output_path folder that stores the output of RaPID.
 *     Assume output of Chromosome i is stored in subfolder i.
 *     E.g. output of Chromosome 1 stored in {rapid_output_path}/1/
 * @param chromosome_start first chromosome to parse.
 * @param chromosome_end last chromosome to parse.
 * @param matrix a 2D unordered map where M[i][j] is a struct pair_stats that
 *     records the total IBD1 and IBD2 between individual with index i and
 *     individual with index j.
 */
static void
worker(
    int thread,
    int num_threads,
    Proceed &proceed,
    Dumpable &dumpable,
    Ordering &order,
    std::string &rapid_output_path,
    int chromosome_start,
    int chromosome_end,
    std::unordered_map<int, std::unordered_map<int, struct pair_stats>> &matrix)
{
    // Input file for each chromosome
    std::vector<std::istream*> inputs;
    chromosome_end = std::min(chromosome_end, NUM_CHROMOSOMES);
    int num_chromosomes = chromosome_end - chromosome_start + 1;
    inputs.reserve(num_chromosomes);

    for (int chrom = chromosome_start; chrom <= chromosome_end; ++chrom) {
        if (chrom > NUM_CHROMOSOMES) {
            continue;
        } else {
            // Read gzipped rapid output file
            std::string file_path = rapid_output_path + "/" + std::to_string(chrom) +"/results.max.gz";
            std::ifstream *file = new std::ifstream(file_path, std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_streambuf<boost::iostreams::input> *buffer = new boost::iostreams::filtering_streambuf<boost::iostreams::input>();
            (*buffer).push(boost::iostreams::gzip_decompressor());
            (*buffer).push(*file);

            inputs.push_back(new std::istream(buffer));
        }
    }

    // Index of the last individual processed for each chromosome
    std::vector<int> prev_ids(num_chromosomes, -1);
    // Segments shared by the last individual for each chromosome
    std::vector<std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<int, int>>>>> chrom_to_id_to_haps_to_segment(num_chromosomes);

    int num_finished_chromosomes = 0;
    std::vector<bool> has_finished(num_chromosomes, false);

    while (num_finished_chromosomes != num_chromosomes) {
        for (int chrom = chromosome_start; chrom <= chromosome_end; ++chrom) {
            int index = chrom - chromosome_start;

            if (has_finished[index]) {
                continue;
            }
            std::istream &in = *inputs[index];
            std::string line;

            int num_finished_ids = prev_ids[index] == -1 ? -1 : 0;
            while (num_finished_ids < NUM_IDS_PER_CYCLE) {
                // Handle NUM_IDS_PER_CYCLE individuals in a cycle
                if (!std::getline(in, line)) {
                    // This chromosome has been exausted
                    has_finished[index] = true;
                    ++num_finished_chromosomes;

                    // Update total IBD1 for the last individual
                    update_total_ibd1(
                        chrom, prev_ids[index],
                        chrom_to_id_to_haps_to_segment[index],
                        matrix
                    );

                    // Discard information about the last individual
                    chrom_to_id_to_haps_to_segment[index].clear();

                    // All individuals can be dumped
                    dumpable.update(chrom, order.get_last_index());

                    break;
                } else {
                    // This chromosome has not been exausted
                    struct line_info info;
                    parse_line(line, info, order);
                    if (info.id1_index == info.id2_index) {
                        continue;
                    }

                    // Segments shared by the individual of the previous segment
                    // on this chromosome
                    std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<int, int>>>> &id_to_haps_to_segment = chrom_to_id_to_haps_to_segment[index];

                    if (process_segment(info, chrom, chromosome_start, prev_ids, id_to_haps_to_segment, matrix)) {
                        // Segments for the previous individual on this chromosome have been exausted
                        ++num_finished_ids;
                        if (num_finished_ids == NUM_IDS_PER_CYCLE) {
                            dumpable.update(chrom, info.id1_index - 1);
                        }
                    }
                }
            }
        }

        {
            std::unique_lock<std::mutex> lock(proceed.mutex);
            if (proceed.increment_num_blocked_and_get() == num_threads) {
                // This thread is the last thread that reaches synchronization.
                // It should signal master thread to dump data.
                proceed.signal_master();
            }

            // Block until master thread has finished dumping data
            while (!proceed.can_thread_proceed(thread)) {
                proceed.wait_for_master(lock);
            }

            // Allow this thread to block at next synchronization
            proceed.disallow_thread_proceed(thread);
        }

    }

    // ALl chromosomes this thread is responsible for have been exausted
    {
        std::unique_lock<std::mutex> lock(proceed.mutex);

        proceed.mark_thread_finished(thread);

        // Last thread to reach synchronization. Signal master thread.
        if (proceed.increment_num_blocked_and_get() == num_threads) {
            proceed.signal_master();
        }
    }
}


/**
 * Update total IBD1 shared between an individual specified by id_index and
 * any other individual using segments on one chromosome.
 *
 * @param chromosome_number
 * @param id_index index of the individual
 * @param id_to_haps_to_segment an unordered map from the index of the other individual
 *     in the pair to a vector of vectors of segments. The each inner vectors contain
 *     segments on haplotypes 0-0, 0-1, 1-0, and 1-1. The exact index is converted using
 *     haps_to_encoding.
 * @param matrix a 2D unordered map where M[i][j]
 *     is a struct pair_stats that records the total IBD1 and IBD2 between individual
 *     with index i and individual with index j.
 */
static void
update_total_ibd1(
    int chromosome_number,
    int id_index,
    std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<int, int>>>> &id_to_haps_to_segment,
    std::unordered_map<int, std::unordered_map<int, struct pair_stats>> &matrix)
{
    for (auto &iter : id_to_haps_to_segment) {
        // Merge segments
        double total_ibd1 = compute_total_ibd1(*merge_four_segment_vectors(
            iter.second[haps_to_encoding(0, 0)],
            iter.second[haps_to_encoding(0, 1)],
            iter.second[haps_to_encoding(1, 0)],
            iter.second[haps_to_encoding(1, 1)]
            ).get(),
            chromosome_number
        );
        // Update total IBD1
        struct pair_stats &pair = matrix[id_index][iter.first];
        pair.total_ibd1 += total_ibd1;
    }
}


static void tokenize_line(const std::string& str,
		std::vector<std::string>& tokens,
		const std::string& delimiters)
{
	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

/**
 * Parse the input line and fill info.
 *
 * @param one line of RaPID output
 * @param info a struct line_info that will be filled
 * @param order an Ordering that specifies the ordering of the IDs as they appear in VCF.
 */
static void
parse_line(
    std::string &line,
    struct line_info &info,
    class Ordering &order)
{
  //  std::vector<std::string> res;//(10);
//    std::cout << "line:" << line << "\n";
    //std::cout << line.find_first_not_of("1", 0);
    std::mutex mtx;           // mutex for critical section
  //mtx.lock();

   // tokenize_line(line,res,"\t");
    std::unique_ptr<char[]> c_line = std::make_unique<char[]>(line.size() + 1);
    strcpy(c_line.get(), line.c_str());
    char *p_save;
   // std::cout << line.c_str();
    char *c_field = strtok_r(c_line.get(), "\t", &p_save);
    c_field = strtok_r(NULL, "\t", &p_save);
    int token_count = 2;
//    std::cout << res[8] << "\n";

  //mtx.unlock();
/**
    info.id1_index = atoi(res[1].c_str());
    info.id2_index = atoi(res[2].c_str());
    info.hap1  = atoi(res[3].c_str());
    info.hap2 = atoi(res[4].c_str());
    info.starting_site = atoi(res[8].c_str());
    info.ending_site =atoi(res[9].c_str());
*/
    //std::cout << info.id1_index << ":id\n";
    //std::cout << info.hap1 << ":hap1\n";
    while (c_field != NULL) {
        std::string field(c_field);

        switch (token_count) {
            case 2:
                // Second field in a line is id1
                info.id1_index = order.get_index(field);
		// std::cout << info.id1_index << ":id\n";
                break;
            case 3:
                // Third field is id2
                info.id2_index = order.get_index(field);
                break;
            case 4:
                // Fourth field is haplotype 1
                info.hap1 = std::stoi(field);
                break;
            case 5:
                // Fifth field is haplotype 2
                info.hap2 = std::stoi(field);
                break;
            case 9:
                // Ninth field is starting site
                info.starting_site = std::stoi(field);
                break;
            case 10:
                // Tenth field is ending site
                info.ending_site = std::stoi(field);
                break;
            default:
                break;
        }


        c_field = strtok_r(NULL, "\t", &p_save);
        ++token_count;
    }

}


/**
 * Handle a newly parsed segment. Update total IBD1 and total IBD2 of the individual.
 *
 * @param info a struct line_info that contains information about the segment.
 * @param chromosome_number the chromosome the segment is on.
 * @param chromosome_start the first chromosome the thread calling this function
 *     is responsible for.
 * @param prev_ids a vector of indices of the last individuals processed on chromosomes
 *     the calling thread is responsible for.
 * @param id_to_haps_to_segment an unordered map from the index of the other individual
 *     in the pair to a vector of vectors of segments. The each inner vectors contain
 *     segments on haplotypes 0-0, 0-1, 1-0, and 1-1. The exact index is converted using
 *     haps_to_encoding.
 * @param matrix a 2D unordered map where M[i][j]
 *     is a struct pair_stats that records the total IBD1 and IBD2 between individual
 *     with index i and individual with index j.
 *
 * @return whether this individual is a new individual on the input chromosome
 */
static bool
process_segment(
    struct line_info &info,
    int chromosome_number,
    int chromosome_start,
    std::vector<int> &prev_ids,
    std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<int, int>>>> &id_to_haps_to_segment,
    std::unordered_map<int, std::unordered_map<int, struct pair_stats>> &matrix)
{
    int hap_encoding = haps_to_encoding(info.hap1, info.hap2);
    int index = chromosome_number - chromosome_start;

    bool is_same_individual_as_last_segment = info.id1_index == prev_ids[index];
    if (!is_same_individual_as_last_segment) {
        // New individual

        update_total_ibd1(chromosome_number, prev_ids[index], id_to_haps_to_segment, matrix);
        id_to_haps_to_segment.clear();

        // Remeber this new individual
        prev_ids[index] = info.id1_index;

    } else {
        // Still the same individual as last one
        const std::unordered_map<int, std::vector<std::pair<int, int>>>::iterator &iter = id_to_haps_to_segment[info.id2_index].find(haps_encoding_to_complement(hap_encoding));
        if (iter != id_to_haps_to_segment[info.id2_index].end()) {
            // Complements exist. Scan them to find IBD2
            std::vector<std::pair<int, int>> &segments = iter->second;
            for (auto &segment : segments) {
                int start = get_intersection_start(info.starting_site, segment.first);
                int end = get_intersection_end(info.ending_site, segment.second);
                if (intersect(start, end)) {
                    // Update total IBD2
                    struct pair_stats &stats = matrix[info.id1_index][info.id2_index];
                    stats.total_ibd2 += get_genetic_length(start, end, chromosome_number);
                }
            }
        }
    }

    // Store this segment.
    id_to_haps_to_segment[info.id2_index][hap_encoding].emplace_back(info.starting_site, info.ending_site);

    return !is_same_individual_as_last_segment;
}


/**
 * @param encoding encoding of a haplotype combination (0-0, 0-1, 1-0, or 1-1).
 *
 * @return encoding of the complement of the combination.
 *     E.g. complement of 0-0 is 1-1.
 */
static inline int
haps_encoding_to_complement(int encoding)
{
    return 3 - encoding;
}


/**
 * @param intersection_start
 * @param intersection_end
 *
 * @return whether the given start and end form a valid segment.
 */
static inline bool
intersect(int intersection_start, int intersection_end)
{
    return intersection_start <= intersection_end;
}


/**
 * @param start1 start of one segment.
 * @param start2 start of another segment.
 *
 * @return start of their intersection assuming they do intersect.
 */
static inline int
get_intersection_start(int start1, int start2)
{
    return std::max(start1, start2);
}


/**
 * @param end1 end of one segment.
 * @param end2 end of another segment.
 *
 * @return end of their intersection assuming they do intersect.
 */
static inline int
get_intersection_end(int end1, int end2)
{
    return std::min(end1, end2);;
}


/**
 * @param hap1 haplotype 1 of a segment.
 * @param hap2 haplotype 2 of the segment.
 *
 * @return an encoding for the haplotype combination hap1-hap2.
 */
static inline int
haps_to_encoding(int hap1, int hap2)
{
    return hap1 + hap2 * 2;
}


/**
 * Merge the four vectors of segments into one vector of segments containing all the old
 * segments sorted by their starts.
 *
 * @param segments1 segments of the first haplotype combination sorted by their starts.
 * @param segments2 segments of the second haplotype combination sorted by their starts.
 * @param segments3 segments of the third haplotype combination sorted by their starts.
 * @param segments4 segments of the fourth haplotype combination sorted by their starts.
 *
 * @return a unique pointer to merged vector of segments.
 */
static inline std::unique_ptr<std::vector<std::pair<int, int>>>
merge_four_segment_vectors(
    std::vector<std::pair<int, int>> &segments1,
    std::vector<std::pair<int, int>> &segments2,
    std::vector<std::pair<int, int>> &segments3,
    std::vector<std::pair<int, int>> &segments4)
{
    return merge_two_segment_vectors(
        *merge_two_segment_vectors(segments1, segments2).get(),
        *merge_two_segment_vectors(segments3, segments4).get()
    );
}

/**
 * Merge the two vectors of segments into one vector of segments containing all the old
 * segments sorted by their starts.
 *
 * @param segments1 segments of the first haplotype combination sorted by their starts.
 * @param segments2 segments of the second haplotype combination sorted by their starts.
 *
 * @return a unique pointer merged vector of segments.
 */
static inline std::unique_ptr<std::vector<std::pair<int, int>>>
merge_two_segment_vectors(
    std::vector<std::pair<int, int>> &segments1,
    std::vector<std::pair<int, int>> &segments2)
{
    int size1 = segments1.size();
    int size2 = segments2.size();
    std::unique_ptr<std::vector<std::pair<int, int>>> p_merged = std::make_unique<std::vector<std::pair<int, int>>>(size1 + size2);
    std::vector<std::pair<int, int>> &merged = *p_merged.get();
    int i = 0;
    int j = 0;
    int fill = 0;
    while (i < size1 && j < size2) {
        int start1 = segments1[i].first;
        int start2 = segments2[j].first;
        if (start1 < start2) {
            merged[fill] = segments1[i];
            ++i;
        } else {
            merged[fill] = segments2[j];
            ++j;
        }
        ++fill;
    }
    while (i < size1) {
        merged[fill] = segments1[i];
        ++i;
        ++fill;
    }

    while (j < size2) {
        merged[fill] = segments2[j];
        ++j;
        ++fill;
    }

    return p_merged;
}


/**
 *
 *
 * @param segments a vector of segments sorted by their starts.
 * @paran chromosome_number the chromosome the segments belong to.
 *
 * @return total IBD1 length of the input segments without overlap.
 */
static inline double
compute_total_ibd1(
    std::vector<std::pair<int, int>> &segments,
    int chromosome_number)
{
    double ans = 0;
    int current_start = segments.front().first;
    int current_end = segments.front().second;
    unsigned int next = 1;

    while (next < segments.size()) {
        int next_start = segments[next].first;
        int next_end = segments[next].second;
        int intersection_start = get_intersection_start(current_start, next_start);
        int intersection_end = get_intersection_end(current_end, next_end);
        if (intersect(intersection_start, intersection_end)) {
            current_end = std::max(current_end, next_end);
        } else {
            ans += get_genetic_length(current_start, current_end, chromosome_number);
            current_start = next_start;
            current_end = next_end;
        }
        ++next;
    }

    ans += get_genetic_length(current_start, current_end, chromosome_number);

    return ans;
}
