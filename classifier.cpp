/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for inferring the relationship between a pair of individuals.
 *
 */

#include <algorithm>
#include <iostream>
#include <unordered_set>

#include <boost/functional/hash.hpp>

#include "parser.hpp"
#include "mapper.hpp"
#include "classifier.hpp"

// Unadjusted thresholds
static double FOURTH_START_ORIGINAL = 1 / std::pow(2, 11.0 / 2);
static double THIRD_START_ORIGINAL = 1 / std::pow(2, 9.0 / 2);
static double SECOND_START_ORIGINAL = 1 / std::pow(2, 7.0 / 2);
static double PO_FS_START_ORIGINAL = 1 / std::pow(2, 5.0 / 2);
static double MZ_START_ORIGNAL = 1 / std::pow(2, 3.0 / 2);
static double FS_START_ORIGINAL = 0.1;

// Current thresholds used for inference
double FOURTH_START = 1 / std::pow(2, 11.0 / 2);
double THIRD_START = 1 / std::pow(2, 9.0 / 2);
double SECOND_START = 1 / std::pow(2, 7.0 / 2);
double PO_FS_START = 1 / std::pow(2, 5.0 / 2);
double MZ_START = 1 / std::pow(2, 3.0 / 2);
double FS_START = 0.1;

// Expected kinship coefficent for full-siblings
#define PO_FS_START_EXPECTED 0.25

const int NUM_TYPES = 7;
const std::vector<std::string> TYPES {"MZ", "PO", "FS", "2nd", "3rd", "4th", "UN"};

// Sum of kinship coefficients of all recorded full-siblings
static double FS_TOTAL_KINSHIP_COEFFICIENTS = 0;
// Number of pairs of full-siblings currently recorded
static int NUM_FS = 0;
// Number of pairs of full-siblings recorded at last adjustment 
static int PREV_ADJUSTED_NUM_FS = 0;

// Mininum number of pairs of full-siblings needed to be recorded between two
// consecutive adjustment.
#define MIN_ADJUSTING_INTERVAL 50

//#define MIN_ADJUSTING_INTERVAL 10

// Maximum number of pairs of full-siblings that will be recorded.
#define MAX_NUM_FS 1000

/**
 * Record an identified full-sibling pair. Add their kinship coefficient to
 * FS_TOTAL_KINSHIP_COEFFICIENTS. If MAX_NUM_FS pairs of
 * full siblings have already been recorded, then no-op.
 *
 * @param stats a struct pair_stats that contains total IBD1 and total IBD2
 *     of the pair.
 *
 */
void
add_full_sibling(const struct pair_stats &stats)
{

// CHANGED! Commented out
    if (NUM_FS > MAX_NUM_FS) {
        return;
    }

    ++NUM_FS;
    FS_TOTAL_KINSHIP_COEFFICIENTS += compute_kinship_coefficient(
        stats.total_ibd1, stats.total_ibd2
    );
}

/**
 *
 * @return number of pairs of full-siblings that have been recorded.
 *
*/
int
get_num_full_siblings() {
    return NUM_FS;
}

/**
 * Adjust kinship coefficient and IBD2 thresholds. No-op if there are fewer than
 * MIN_NUM_FS pairs of full-siblings, or fewer than MIN_ADJUSTING_INTERVAL pairs
 * of full-siblings have been recorded since last adjustment, or there more more
 * than MAX_NUM_FS pairs of full-siblings.
 *
 */
void
shift_boundary()
{

    if (NUM_FS - PREV_ADJUSTED_NUM_FS < MIN_ADJUSTING_INTERVAL ||
        NUM_FS < MIN_NUM_FS || NUM_FS > MAX_NUM_FS) {
        return;
    }


/**
    if (NUM_FS - PREV_ADJUSTED_NUM_FS < MIN_ADJUSTING_INTERVAL ||
        NUM_FS < MIN_NUM_FS) {
        return;
    }
**/


    double mean = FS_TOTAL_KINSHIP_COEFFICIENTS / NUM_FS;
    double shift = mean / PO_FS_START_EXPECTED;

    if ( shift > 1){
	shift = 1;
	}
    std::cout << "Shifting factor is " << shift << std::endl;


    FOURTH_START = FOURTH_START_ORIGINAL * shift;
    THIRD_START = THIRD_START_ORIGINAL * shift;
    SECOND_START = SECOND_START_ORIGINAL * shift;
    PO_FS_START = PO_FS_START_ORIGINAL * shift;
    MZ_START = MZ_START_ORIGNAL * shift;
    FS_START = FS_START_ORIGINAL * shift;

    PREV_ADJUSTED_NUM_FS = NUM_FS;
}


