#ifndef _JEWELER_H_
#define _JEWELER_H_

#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <map>
#include <set>
#include <vector> 
#include "fasta.hpp"
#include "common.hpp"

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include <Fasta.h>

using namespace BamTools;


#include <boost/filesystem.hpp>

using boost::filesystem::path;
using boost::filesystem::create_directory;

class SewingMachine;
class JewelerInfo;
class jeweler{
public:
	
	string info_filename;
	string mamf_filename;
	string mismatch_filename;
	string bracelet_filename;

	JewelerInfo * jeweler_info ;
	FILE * log_file;

	SewingMachine *sm;
	jeweler(int argc, char *argv[]);

	// debug
	// test case 
	int test_case;
    int num_test_case;

	bool is_earrings;
	bool is_bracelet;
	bool is_mismatch_analyzer;
	bool is_prepare;
	// function to run the analysis 
	int run();

	void load_mamf_file();
	
};

#endif /* _JEWELER_H_ */
