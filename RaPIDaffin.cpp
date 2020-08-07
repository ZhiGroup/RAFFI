/**
 * Author: Junjie Shi, Ardalan Naseri
 *
 * This file contains entry point of the program.
 *
 */

#include "RaPIDaffin.hpp"
#ifndef NULL
#define NULL   ((void *) 0)
#endif
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <memory>
#include <string>
#include <boost/thread.hpp>
#include <cmath>
#include <unordered_set>
#include <vector>
#include <boost/filesystem.hpp>
#include <sys/wait.h>
#include "mapper.hpp"
#include "parser.hpp"
#include "classifier.hpp"
#include <sstream>
using namespace std;

#define NUM_REQUIRED_OPTIONS 3

static bool parse_parameters(int args, char** argv, struct parameter &parameters);
static void print_usage(std::ostream &out);

struct parameter {
	std::string input_folder_vcf_path;
	std::string vcf_prefix;
	std::string gen_map_path;
	std::string output_path;
	std::string vcf_example;
	std::string python_path="python3.6";

	std::string rapid_output_path;
	int rapid_out_put_set = 0;
	int max_degree = 4;
	unsigned int num_threads = 22;
};



static std::string exec(const char* cmd) {
	std::array<char, 128> buffer;
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
	if (!pipe) {
		throw std::runtime_error("popen() failed!");
	}
	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
		result += buffer.data();
	}
	return result;
}

int main(int args, char** argv)
{
	struct parameter params;
	if (!parse_parameters(args, argv, params)) {
		print_usage(std::cerr);
		return -1;
	}

	if (!params.output_path.empty()) {
		try {
			boost::filesystem::create_directories(params.output_path);
		} catch (boost::filesystem::filesystem_error &e) {
			std::cerr << "Failed to create output directory: \"" + params.output_path + "\"!" << std::endl;
			return -1;
		}
	}

	stringstream input_vcf_file_example;
	input_vcf_file_example << params.input_folder_vcf_path + "//" + params.vcf_prefix  << 22 <<".vcf.gz";
	string iv = input_vcf_file_example.str();
	params.vcf_example = iv.c_str();

	//cout << params.vcf_example << "\n";
	if (params.rapid_out_put_set == 0){

	//Run RaPID
	stringstream command_line;
	command_line  <<params.python_path<< " " << "../bin/estimate_params.py ";
	command_line << params.vcf_example << " " << params.gen_map_path << "/" << "chr22.rMap";

	//cout << command_line.str() <<"\n";
	int window_size = stoi(exec(command_line.str().c_str()));

	int chr_per_thread = 22 / params.num_threads;
	if (chr_per_thread < 1) chr_per_thread = 1;
	stringstream rapid_params;
	rapid_params << "../bin/RaPID_v.1.7 -r 3 -s 1 -d 5 ";
	rapid_params << " -w " << window_size;
	int chr_counter = 1;
	vector<pid_t> pids(params.num_threads);

	vector<string> clines;//(params.num_threads);

	while (chr_counter <= 22){
	stringstream rapid_command_line;

		for (int j = 0; j < chr_per_thread ; j++){
			if (chr_counter == 23) break;

			rapid_command_line << rapid_params.str() << " -i " << params.input_folder_vcf_path <<"/" <<params.vcf_prefix << chr_counter << ".vcf.gz ";
			rapid_command_line << " -o " << params.output_path << "/" << chr_counter;
			rapid_command_line << " -g " << params.gen_map_path << "/" << "chr" << chr_counter<<".rMap";

			if (j < chr_per_thread -1 )
				rapid_command_line << ";";
			chr_counter++;

		}
		clines.push_back(rapid_command_line.str());
	}


	for(unsigned int i=0;i<params.num_threads;i++) // loop will run n times (n=5)
	    {		        	//cout << clines[i].c_str() <<"\n";


	        if(fork() == 0)
	        {
	        	cout << clines[i].c_str() <<"\n";
	        	exec(clines[i].c_str());
	        	exit(0);
	        }
	    }
	    for(unsigned int i=0;i<params.num_threads;i++)
	    wait(NULL);
	    params.rapid_output_path = params.output_path;

	}

	std::ofstream out(params.output_path + "predictions.txt", std::ios::out | std::ios::trunc);
	out << std::fixed << std::setprecision(4);


	{
	master(
			params.vcf_example, params.rapid_output_path, params.gen_map_path,
			params.max_degree, params.num_threads, out
	);
	}

	return 0;
}


static void
print_usage_old(std::ostream &out)
{
	out << "Usage: ./RaAFFI_v.0.1" << std::endl
			<< "Required parameters:" << std::endl
			<< "-v {vcf path}" << std::endl
			<< "\tPath to any of the gzipped VCF files that were given to RaPID as inputs." << std::endl
			<< "-i {rapid output path} " << std::endl
			<< "\tDirectory containing RaPID output directories." << std::endl
			<< "\tOutput for Chromosome i must be in subdirectory {i}." << std::endl
			<< "-m {map path}" << std::endl
			<< "\tDirectory containing genetic maps that were given to RaPID as inputs." << std::endl
			<< "\tMap for Chromosome i must be named chr{i}.rMap." << std::endl
			<< "Optional parameters:" << std::endl
			<< "-o {output directory}" << std::endl
			<< "\tOutput will be written to {output directory}/predictions.txt." << std::endl
			<< "\tDefault is current direcotry." << std::endl
			<< "-d {max degree}" << std::endl
			<< "\tMaximum target degree (4 is largest supported degree)." << std::endl
			<< "\tDefault is 4." << std::endl
			<< "-p {python version}" << std::endl
			<< "\tOptimal number of threads is the number of chromosomes." << std::endl
			<< "\tDefault is 22." << std::endl;


}

static void
print_usage(std::ostream &out)
{
	out << "Usage: ./RAFFI_v.0.1" << std::endl
			<< "Required parameters:" << std::endl
			<< "-i <input folder>" << std::endl
			<< "\tPath to the folder containing gzipped VCF files.." << std::endl
			<< "-v <vcf prefix> " << std::endl
			<< "\tPrefix of gzipped VCF files. e.g. hap1_chr for hap1_chr1.vcf.gz" << std::endl
			<< "-g {genetic mapping file path}" << std::endl
			<< "\tDirectory containing genetic maps that were given to RaPID as inputs." << std::endl
			<< "\tMap for Chromosome i must be named chr{i}.rMap." << std::endl

			<< "Optional parameters:" << std::endl
			<< "-O {output directory of RaPID results}" << std::endl
			<< "\tUsing this tag, RaPID will not be run and the provided outputs will be used directly for relatedness inference." << std::endl
			<< "-o {output directory}" << std::endl
			<< "\tOutput will be written to {output directory}/predictions.txt." << std::endl
			<< "\tDefault is current direcotry." << std::endl
			<< "-d {max degree}" << std::endl
			<< "\tMaximum target degree (4 is largest supported degree)." << std::endl
			<< "\tDefault is 4." << std::endl
			<< "-t {number of threads}" << std::endl
			<< "\tOptimal number of threads is the number of chromosomes." << std::endl
			<< "\tDefault is 22." << std::endl
		    << "-p {Python version}" << std::endl
		    << "\tPython path" << std::endl
			<< "\tDefault is python3.6" << std::endl;
}

/**
 * Parse command line arguments and fill parameters.
 *
 * @param args number of command line arguments.
 * @param array of command line arguments.
 * @param parameters a struct parameter that will be filled.
 *
 * @return whether the arguments are valid and required parameters not missing.
 */
static bool
parse_parameters(int args, char** argv, struct parameter &parameters)
{
	int i = 1;
	bool failed = false;
	std::unordered_set<std::string> detected_options;
	while (i < args) {
		std::string option = argv[i];
		if (option == "-v") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			detected_options.insert(option);
			parameters.vcf_prefix = argv[i];
		//	parameters.vcf_prefix += "/";
		} else if (option == "-O") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			parameters.rapid_out_put_set = 1;
			parameters.rapid_output_path = argv[i];
			parameters.rapid_output_path += "/";
		}else if (option == "-o") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			parameters.output_path = argv[i];
			parameters.output_path += "/";
		}

		else if (option == "-i") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			detected_options.insert(option);
			parameters.input_folder_vcf_path = argv[i];
		} else if (option == "-g") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			detected_options.insert(option);
			parameters.gen_map_path = argv[i];
			parameters.gen_map_path += "/";
		} else if (option == "-d") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			parameters.max_degree = std::stoi(argv[i]);
			if (parameters.max_degree > 4 || parameters.max_degree <= 0) {
				std::cerr << "Degrees less than 1 or beyond 4 are not supported!" << std::endl << std::endl;
				failed = true;
				break;
			}
		} else if (option == "-t") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			parameters.num_threads = std::stoul(argv[i], nullptr);
		}else if (option == "-p") {
			i++;
			if (i >= args) {
				failed = true;
				break;
			}
			detected_options.insert(option);
			parameters.python_path = argv[i];
		}
		i++;
	}
	if (detected_options.size() < NUM_REQUIRED_OPTIONS) {
		failed = true;
	}
	return !failed;
}
