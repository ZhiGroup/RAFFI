/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for genetic map related tasks.
 *
 */

#ifndef MAPPER_HPP
#define MAPPER_HPP

#include <string>

extern double TOTAL_LENGTH;

void init_maps(std::string &map_path);
double get_genetic_length(int starting_site, int ending_site, int chromosome_number);
void deinit_maps();

#endif
