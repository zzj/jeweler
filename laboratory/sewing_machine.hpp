#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <map>
#include <Fasta.h>
using namespace BamTools;
using namespace std;

#include "../common.hpp"
#include "bam_info.hpp"
#include "cigar_holder.hpp"
class locator{
public:
	int ref_id;
	int position;
	int num_mismatches;
	int edit_distance;
	bool is_first_pair;
	string cigar_string;
	locator(BamAlignment &al);
};

typedef vector<locator *> locator_pvector;


class compare_locator {
public:
	bool operator()(const locator* a, const locator *b);
};

class SewingMachine: public BamInfo{
private:

	int output_alignment_map_core(FILE *, bool only_multiple_reads);
	int output_alignment_connection_map_core(FILE *);
	int build_alignment_connection_map_core();

public:
	
	SewingMachine();
	SewingMachine(BamReader &reader);
	int initialize(BamReader &reader);
	
	// functions
	int add_alignment(BamAlignment &al);
	int output_alignment_map(FILE *);
	int output_multiple_alignment_map(FILE *);
	int output_alignment_connection_map(FILE *);
	string get_reference_sequence(FastaReference &fr, BamAlignment &al);
	// storage
	map<string,locator_pvector> seqs;

	// statistics 
	vector<int> num_multiple_reads_per_chr;
	
	// map
	map<locator *, locator_pvector, compare_locator> alignment_connection_map;
	
};
