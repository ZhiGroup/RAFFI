/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for parsing the output of RaPID.
 *
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include <unordered_map>
#include <memory>

void master(
    std::string &vcf_path,
    std::string &rapid_output_path,
    std::string &map_path,
    int max_degree,
    unsigned int num_threads,
    std::ostream &out);

struct pair_stats {
        double total_ibd1 = 0;
        double total_ibd2 = 0;
};

#endif
