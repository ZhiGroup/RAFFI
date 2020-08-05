/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 *
 *
 */

#ifndef DUMPABLE_HPP
#define DUMPABLE_HPP

#include <vector>
#include <iostream>

#include "ordering.hpp"
#include "classifier.hpp"


void infer_candidates(
    double max_degree,
    uint64_t num_dumped, std::istream &temp_in,
    Ordering &order, std::ostream &out);


inline void
output_header(std::ostream &out)
{
    out << "ID1\tID2\tKINSHIP\tIBD0\tIBD1\tIBD2\tTYPE" << std::endl;
}


inline void
write_pair(
    std::string &id1, std::string &id2,
    double kinship_coefficient, double probability_ibd0,
    double probability_ibd1, double probability_ibd2,
    int encoding,
    std::ostream &out)
{
    out << id1 << "\t" << id2 << "\t"
        << kinship_coefficient << "\t"
        << probability_ibd0 << "\t"
        << probability_ibd1 << "\t"
        << probability_ibd2 << "\t"
        << TYPES[encoding] << std::endl;
}


int dump_range(
    int max_degree,
    double min_kinship_coefficient,
    Ordering &order,
    const std::pair<int, int> &range,
    std::vector<std::unordered_map<int, std::unordered_map<int, struct pair_stats>>> &matrices,
    std::ostream &temp_out,
    std::ostream &out);


struct dumpable_pair {
    int id1_index;
    int id2_index;
    double kinship_coefficient;
    double probability_ibd2;
};


class Dumpable {
public:
    Dumpable(const class Ordering &order);

    void update(int chromosome_number, int index);

    std::pair<int, int> get_dumpable_indices();

    void set_previous_last_dumpable_index(int index);

private:
    int previous_last_dumpable_index;
    std::vector<int> last_dumpable_indices;
    const class Ordering &id_ordering;
};



#endif

