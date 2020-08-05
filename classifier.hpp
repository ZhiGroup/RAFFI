/**
 * Author: Junjie Shi (2019)
 *
 * This file is responsible for infering the relationship between a pair of individuals.
 *
 */

#ifndef CLASSIFIER_HPP
#define CLASSIFIER_HPP

#include <vector>
#include <string>
#include <unordered_map>

#include "mapper.hpp"

const extern int NUM_TYPES;
const extern std::vector<std::string> TYPES;

extern double FOURTH_START;
extern double THIRD_START;
extern double SECOND_START;
extern double PO_FS_START;
extern double MZ_START;
extern double FS_START;

#define MIN_NUM_FS 200

void add_full_sibling(const struct pair_stats &stats);
void shift_boundary();
int get_num_full_siblings();


inline int
get_encoding(double kinship_coefficient, double probability_ibd2)
{
   if (FOURTH_START <= kinship_coefficient && kinship_coefficient < THIRD_START) {
      return 5;
   } else if (THIRD_START <= kinship_coefficient && kinship_coefficient <= SECOND_START) {
      return 4;
   } else if (SECOND_START <= kinship_coefficient && kinship_coefficient < PO_FS_START) {
      return 3;
   } else if (PO_FS_START <= kinship_coefficient && kinship_coefficient < MZ_START) {
      return probability_ibd2 >= FS_START ? 2 : 1;
   } else if (MZ_START <= kinship_coefficient) {
      return 0;
   } else {
      return NUM_TYPES - 1;
   }
}


inline double
compute_kinship_coefficient_from(
   double probability_ibd1, double probability_ibd2)
{
   return probability_ibd1 / 4 + probability_ibd2 / 2;
}

inline double
compute_kinship_coefficient(
   double total_ibd1, double total_ibd2)
{
   return total_ibd1 / (4 * TOTAL_LENGTH) + total_ibd2 / (2 * TOTAL_LENGTH);
}


inline double
compute_probability_ibd0(double total_ibd1, double total_ibd2)
{
   return (TOTAL_LENGTH - total_ibd1 - total_ibd2) / TOTAL_LENGTH;
}

inline double
compute_probability_ibd1(double total_ibd1)
{
   return total_ibd1 / TOTAL_LENGTH;
}

inline double
compute_probability_ibd2(double total_ibd2)
{
   return total_ibd2 / TOTAL_LENGTH;
}


#endif
