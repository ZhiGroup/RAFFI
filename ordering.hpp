/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for determining the ordering of the IDs in the VCF files.
 *
 */

#ifndef ORDERING_HPP
#define ORDERING_HPP

#include <unordered_map>
#include <vector>
#include <string>




class Ordering {
	std::unordered_map<std::string, int> id_to_index;
	std::vector<std::string> ids;

public:
	Ordering(std::string &file_path);

	int get_index(std::string &id);

	int get_last_index();

	std::string& get(int index);

	int size();



private:
	void get_id_ordering(std::string &file_path);
};

#endif
