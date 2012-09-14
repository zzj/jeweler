#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <Fasta.h>

using namespace BamTools;
using namespace std;

class Appraiser {
public:
	Appraiser(int argc, char *argv[]);

	string bamfile;
	string fastafile; 
	FILE *log_file; // log file
	string mam_table_file; // multiple alignment file
	string mam_map_file; // multiple alignment file
	FILE *quality_file; //quality file
	int test_num;

	void run();

};


