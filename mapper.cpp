/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for genetic map related tasks.
 *
 */

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include "mapper.hpp"

#include "RaPIDaffin.hpp"


// Length of genome
double TOTAL_LENGTH = 0;
// Genetic maps for all chromosomes
static std::vector<std::unique_ptr<std::vector<double>>> maps;

static std::unique_ptr<std::vector<double>> parse_map(int chromosome_number, std::string &map_path);


/**
 * Read genetic maps for all the chromosomes.
 *
 * @param map_path folder in which genetic maps are stored.
 */
void
init_maps(std::string &map_path)
{
    for (int chrom = 1; chrom <= NUM_CHROMOSOMES; ++chrom) {
        maps.push_back(parse_map(chrom, map_path));
    }
}


/**
 * Free all read genetic maps.
 */
void
deinit_maps()
{
    maps.clear();
}


template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * Read the genetic map of one chromosome. Add the length of the chromosome to
 * TOTAL_LENGTH.
 *
 * @param chromosome_number
 * @param map_path folder in which the genetic map is stored. Assume the the map
 *     is named as chr{xx}.rMap.
 *
 * @return unique pointer to a vector of genetic distances.
 */
static std::unique_ptr<std::vector<double>>
parse_map(int chromosome_number, std::string &map_path)
{

    //std::cout << chromosome_number << "," << map_path << "\n";
    std::unique_ptr<std::vector<double>> p_distances = make_unique<std::vector<double>>();
    std::vector<double> &distances = *p_distances.get();

    std::ifstream in {map_path + "chr" + std::to_string(chromosome_number) + ".rMap"};
    std::string line;

    while (std::getline(in, line)) {
        int sep = line.find_first_of("\t");

        double distance = std::strtod(line.substr(sep + 1).c_str(), NULL);
        distances.push_back(distance);
    }

    TOTAL_LENGTH += distances.back() - distances.front();
    return p_distances;
}


/**
 * Compute the genetic length between two sites.
 *
 * @param starting_sitesite
 * @param ending_site
 *
 * @return genetic length
 */
double
get_genetic_length(int starting_site, int ending_site, int chromosome_number)
{
    std::vector<double> &distances = *maps[chromosome_number - 1].get();
    return distances[ending_site] - distances[starting_site];
}

