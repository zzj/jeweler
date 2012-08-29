
#ifndef _SEWING_MACHINE_
#define _SEWING_MACHINE_

#include <iostream>
#include <string>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <map>
#include <set>
#include <Fasta.h>
using namespace BamTools;
using namespace std;

#include "../common.hpp"
#include "bam_info.hpp"
#include "cigar_holder.hpp"

#include "gtest/gtest_prod.h"

class locator{
public:
	int ref_id;
	int position;
	int num_mismatches;
	int edit_distance;
	bool is_first_pair;
	string cigar_string;
	locator(const JewelerAlignment &al);
};

typedef vector<locator *> locator_pvector;


class compare_locator {
public:
	bool operator()(const locator* a, const locator *b);
};

/* this class is used for generating the mulitple alignment map */

class SewingMachine: public BamInfo{
private:

    friend class SewingMachineTest;
    FRIEND_TEST(SewingMachineTest, test_get_full_name);
    FRIEND_TEST(SewingMachineTest, test_count_alignment);
    FRIEND_TEST(SewingMachineTest, test_add_alignment);

	int output_alignment_map_core(FILE *, bool only_multiple_reads);
	int output_alignment_connection_map_core(FILE *);
	int build_alignment_connection_map_core();
    string get_full_name(const JewelerAlignment &al);
public:

	SewingMachine();
	SewingMachine(BamReader &reader);
	void initialize(BamReader &reader);
	void load_multiple_alignments_set(FILE *);
	// functions
	// first time count_alignment
	// second time add_alignment
	int count_alignment(const JewelerAlignment &al);
    int add_alignment(const JewelerAlignment &al);
	int output_alignment_map(FILE *);
	int output_multiple_alignment_map(FILE *);
	int output_alignment_connection_map(FILE *);
	string get_reference_sequence(FastaReference &fr, JewelerAlignment &al);
	// storage
	map<string,locator_pvector> seqs;

	map<string, int> multiple_alignment_map;
	set<string> multiple_alignment_set;

	// statistics
	vector<int> num_multiple_reads_per_chr;

	// map
	map<locator *, locator_pvector, compare_locator> alignment_connection_map;

};

#endif
