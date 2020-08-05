/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file is responsible for determining the ordering of the IDs in the VCF files.
 *
 */

#include <fstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "ordering.hpp"

#include <iostream>

/**
 * Constructor of Ordering.
 *
 * @param file_path path to an VCF file.
 */
Ordering::Ordering(std::string &file_path)
{
    get_id_ordering(file_path);
}

/**
 * @param id an ID
 *
 * @return the index of the ID in the ordering
 */
int
Ordering::get_index(std::string &id)
{
    return id_to_index[id];
}


/**
 * @return index of the last invidual.
 */
int
Ordering::get_last_index()
{
    return ids.size() - 1;
}


/**
 * @index index of an individual
 *
 * @return ID of the individual
 */
std::string&
Ordering::get(int index)
{
    return ids[index];
}


/**
 * @return total number of individuals.
 */
int
Ordering::size()
{
    return id_to_index.size();
}


/**
 * Determines the ordering of the IDs in an VCF file.
 *
 * @param vcf_path path to the VCF file.
 */
void
Ordering::get_id_ordering(std::string &vcf_path)
{
    std::ifstream file(vcf_path, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> buffer;
    buffer.push(boost::iostreams::gzip_decompressor());
    buffer.push(file);
    std::istream in(&buffer);

    std::string line;
    int index = 0;
    while (std::getline(in, line)) {
        if (line.rfind("##", 0) == 0) {
            continue;
        }

        std::unique_ptr<char[]> c_line = std::make_unique<char[]>(line.size() + 1);
        strcpy(c_line.get(), line.c_str());
        char *p_save;
        char *field = strtok_r(c_line.get(), "\t", &p_save);
        int token_count = 1;

        while (field != NULL) {
            // IDs start from 10th field
            if (token_count > 9) {
                id_to_index[field] = index;
                ids.push_back(std::move(field));
                index++;
            }

            field = strtok_r(NULL, "\t", &p_save);
            token_count++;
        }
        // Only need to read the header line
        break;
    }
}
