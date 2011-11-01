#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <map>
using namespace BamTools;
using namespace std;

#include "../common.hpp"
#include "bam_info.hpp"
#include "cigar_holder.hpp"
class locator{
public:
	int ref_id;
	int position;
	bool is_first_pair;
	string cigar_string;
	locator(int ref_id, int position, string cigar_string);
};

typedef vector<locator> locator_vector;

class SewingMachine: public BamInfo{
public:
	SewingMachine();
	SewingMachine(BamReader &reader);


	// functions
	int add_alignment(BamAlignment &al);
	int output_alignment_map(FILE *);
	int output_multiple_alignment_map(FILE *);
	// storage
	map<string,locator_vector> seqs;
	
};
