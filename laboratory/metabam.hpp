
#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "../common.hpp"
#include "bam_info.hpp"
using namespace BamTools;
using namespace std;

class Metabam: public BamInfo{
public:
	Metabam();
	Metabam(BamReader &reader);
	int initialize(BamReader &reader);

	// functions
	int add_alignment(BamAlignment &al);
	int dump_meta_data(FILE * log_file,bool is_readable=true);
	int dump_mapquality_list(FILE *output_file);
	// statistics 
	int total_num_alignments;
	int total_num_duplicate;
	int total_num_first_mate;
	int total_num_second_mate;
	int total_num_mapped;
	int total_num_mate_mapped;
	int total_num_paired;
	int total_num_proper_pair;
	int total_num_reverse_strand;
	vector<int> num_alignment_per_reference;
	vector<int> mapquality_list;
};
