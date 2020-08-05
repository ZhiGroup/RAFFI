/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 *
 * This class is used to write the output results of relatedness inference based on the given max degree.
 *
 */

#include <algorithm>
#include <unordered_map>

#include <boost/iostreams/filtering_streambuf.hpp>

#include "parser.hpp"
#include "dumpable.hpp"
#include "RaPIDaffin.hpp"

static inline bool is_encoding_less_than(int encoding, int degree);
static inline double compute_probability_ibd1_from(
    double kinship_coefficient, double probability_ibd2);

/**
 * Constructor of Dumpable.
 *
 * @param order an Ordering that specifies the ordering of the IDs as they appear
 *     in the vcf file.
 */
Dumpable::Dumpable(const class Ordering &order) :
    previous_last_dumpable_index(-1),
    last_dumpable_indices(NUM_CHROMOSOMES, -1),
    id_ordering(order) {}


/**
 * Notify this Dumpable the index of the last individual as specified by the Ordering
 * that has been processed for the input chromosome. In other words, information about
 * this individual and all individuals before him is no longer needed for processing
 * remaining individuals for this chromosome.
 *
 * @param chromosome_number which chromosome.
 * @param index index of the individual
 */
void
Dumpable::update(int chromosome_number, int index)
{
    last_dumpable_indices[chromosome_number - 1] = index;
}


/**
 *
 * @return inclusive range of indices where all individuals in it have been processed
 *     across all chromosomes and therefore can be written to output.
 *
 */
std::pair<int, int>
Dumpable::get_dumpable_indices()
{
    return {
        previous_last_dumpable_index + 1,
        *std::min_element(last_dumpable_indices.begin(), last_dumpable_indices.end())
    };
}


/**
 * Next call to Dumpable::get_dumpable_indices is guaranteed to return a range
 * starting from the input index plus one.
 *
 * @param the last index in the range returned by a previous call to
 *     Dumpable::get_dumpable_indices.
 */
void
Dumpable::set_previous_last_dumpable_index(int index)
{
    previous_last_dumpable_index = index;
}


/**
 * Write all individuals in the given inclusive range to either temporary output or
 * final output. An individual will be written to temporary output if there are
 * insufficient number of pairs of full-siblings recorded.
 *
 * @param max_degree largest degree user is looking for.
 * @min_kinship_coefficient pairs with kinship coefficients below this will not be
 *     written to any output.
 * @order an Ordering that specifies the ordering of the IDs as they appear in VCF.
 * @range inclusive range of indices in which all individuals can be written.
 * @matrices an vector of matrix. Each matrix is a 2D unordered map where M[i][j]
 *     is a struct pair_stats that records the total IBD1 and IBD2 between individual
 *     with index i and individual with index j.
 * @temp_out temporary output.
 * @out final output.
 *
 * @return number of individuals written to output.
 */
int
dump_range(
    int max_degree,
    double min_kinship_coefficient,
    Ordering &order,
    const std::pair<int, int> &range,
    std::vector<std::unordered_map<int, std::unordered_map<int, struct pair_stats>>> &matrices,
    std::ostream &temp_out,
    std::ostream &out)
{
    int num_dumped = 0;

    for (int id1_index = range.first; id1_index <= range.second; ++id1_index) {
        std::unordered_map<int, struct pair_stats> id2_index_to_writable_stats;

        // Aggregate individuals sharing IBD with this individual
        for (const auto &matrix : matrices) {
            const auto &id1_index_to_id2_index_to_stats = matrix.find(id1_index);
            if (id1_index_to_id2_index_to_stats != matrix.end()) {
                for (const auto &id2_index_to_stats : id1_index_to_id2_index_to_stats->second) {

                    int id2_index = id2_index_to_stats.first;
                    const struct pair_stats &sub_stats = id2_index_to_stats.second;

                    struct pair_stats &stats = id2_index_to_writable_stats[id2_index];
                    stats.total_ibd1 += sub_stats.total_ibd1 - sub_stats.total_ibd2;
                    stats.total_ibd2 += sub_stats.total_ibd2;
                }
            }
        }

        // Write out all the pairs consisted of this individual and another individual
        // sharing IBD with this individual
        for (const auto &id2_index_to_stats : id2_index_to_writable_stats) {
            const struct pair_stats &stats = id2_index_to_stats.second;
            double kinship_coefficient = compute_kinship_coefficient(stats.total_ibd1, stats.total_ibd2);
            double probability_ibd2 = compute_probability_ibd2(stats.total_ibd2);

            if (probability_ibd2 >= FS_START) {
                add_full_sibling(stats);
            }

            int num_full_siblings = get_num_full_siblings();
            if (kinship_coefficient >= min_kinship_coefficient && num_full_siblings < MIN_NUM_FS) {
                // Write to temporary
                struct dumpable_pair pair;
                pair.id1_index = id1_index;
                pair.id2_index = id2_index_to_stats.first;
                pair.kinship_coefficient = kinship_coefficient;
                pair.probability_ibd2 = probability_ibd2;

                if (!temp_out.write(reinterpret_cast<char *>(&pair), sizeof(struct dumpable_pair))) {
                    throw std::runtime_error {"Failed to write to temporary file"};
                }
                ++num_dumped;

            } else if (num_full_siblings >= MIN_NUM_FS) {
                // Write to final directly

                // Adjust inference boundaries
                shift_boundary();

                int encoding = get_encoding(kinship_coefficient, probability_ibd2);
                if (is_encoding_less_than(encoding, max_degree)) {
                    double probability_ibd1 = compute_probability_ibd1(stats.total_ibd1);
                    double probability_ibd0 = std::max(1 - probability_ibd1 - probability_ibd2, 0.0);

                    write_pair(
                        order.get(id1_index),
                        order.get(id2_index_to_stats.first),
                        kinship_coefficient,
                        probability_ibd0,
                        probability_ibd1,
                        probability_ibd2,
                        get_encoding(kinship_coefficient, probability_ibd2),
                        out
                    );
                }
            }
        }

        // Discard information about this individual
        for (auto &matrix : matrices) {
            matrix.erase(id1_index);
        }
    }

    out.flush();
    temp_out.flush();

    return num_dumped;
}


/**
 * Read from temporary output, infer the relationships of the candidates, and write
 * them to final output.
 *
 * @param largest degree user is looking for.
 * @num_dumped total number of pairs written to temporary output.
 * @temp_in temporary output for reading.
 * @order an Ordering that specifies the ordering of the IDs as they appear in VCF.
 * @out final output.
 */
void
infer_candidates(
    double max_degree,
    uint64_t num_dumped, std::istream &temp_in,
    Ordering &order, std::ostream &out)
{
    struct dumpable_pair *pair;
    char buffer[sizeof(struct dumpable_pair)];
    uint64_t num_read = 0;

    // Read a struct dumpable_pair into buffer
    while (temp_in.read(buffer, sizeof(struct dumpable_pair))) {
        ++num_read;
        pair = reinterpret_cast<struct dumpable_pair*>(buffer);

        int encoding = get_encoding((*pair).kinship_coefficient, (*pair).probability_ibd2);

        // Write to final output if this candidate pair is closer than max_degree
        if (is_encoding_less_than(encoding, max_degree)) {
            double probability_ibd1 = std::max(compute_probability_ibd1_from((*pair).kinship_coefficient, (*pair).probability_ibd2), 0.0);
            double probability_ibd0 = std::max(1 - probability_ibd1 - (*pair).probability_ibd2, 0.0);

            write_pair(
                order.get((*pair).id1_index),
                order.get((*pair).id2_index),
                (*pair).kinship_coefficient,
                probability_ibd0,
                probability_ibd1,
                (*pair).probability_ibd2,
                encoding,
                out
            );
        }
    }

    if (num_read != num_dumped || !temp_in.eof()) {
        throw std::runtime_error {"Failed to read from temporary file"};
    }

    out.flush();
}

/**
 * Determines if the input encoding represents a degree closer than the input degree.
 *
 * @param encoding encoding of the relationship between a pair of individuals.
 * @param degree the degree to compare against.
 */
static inline bool
is_encoding_less_than(int encoding, int degree)
{
    if (degree == 1) {
        return encoding <= 2;
    } else {
        return encoding - 1 <= degree;
    }
}

/**
 * Compute probability of IBD1 from kinship coefficient and probability of IBD2.
 *
 * @param kinship_coefficient
 * @param probability_ibd2
 *
 * @return probability of IBD1
 */
static inline double
compute_probability_ibd1_from(double kinship_coefficient, double probability_ibd2)
{
    return (4 * kinship_coefficient * TOTAL_LENGTH - 2 * probability_ibd2) / TOTAL_LENGTH;
}
