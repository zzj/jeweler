
#ifndef _SEWING_MACHINE_
#define _SEWING_MACHINE_

#include <iostream>
#include <string>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <map>
#include <set>
#include <Fasta.h>
#include <memory>
using namespace BamTools;
using namespace std;

#include "../proto/jeweler.pb.h"
#include "../common.hpp"
#include "bam_info.hpp"
#include "cigar_holder.hpp"
#include "gtest/gtest_prod.h"

class ZLevelDB;
class locator{
public:
	locator(const JewelerAlignment &al, const RefVector &rv);
    void dump(Jeweler::SewingMachineData::Locator *);
    bool is_first();
private:
    friend class compare_locator;
    Jeweler::SewingMachineData::Locator data;
};

class compare_locator {
public:
	bool operator()(const locator* a, const locator *b);
};

/* this class is used for generating the mulitple alignment map */

class SewingMachine: public BamInfo{
private:

    friend class SewingMachineTest;

public:

	SewingMachine();
    ~SewingMachine();	
	SewingMachine(BamReader &reader, string &db_file);
	void initialize(BamReader &reader, string &db_file);
    int add_alignment(const JewelerAlignment &al);
    void add_read(const string &name,
                  shared_ptr<Jeweler::SewingMachineData> smd);
    void load_zdb(string &db_file);
	string get_reference_sequence(FastaReference &fr, JewelerAlignment &al);
    bool is_multiple_alignment(const shared_ptr<Jeweler::SewingMachineData> smd);
    bool is_multiple_alignment(const string &name) const;
    int build_alignment_connection_map_core();
    int output_alignment_connection_map_core(FILE * file);
    int output_alignment_connection_map(FILE * file);
    void snapshot();
    string snapshot_file();
    void load_snapshot();
private:
	// storage
    ZLevelDB *zdb;
    string db_folder;
    set<string> multiple_alignment_names;
    size_t initialized_size;
	// statistics
	//vector<int> num_multiple_reads_per_chr;

	// map
	//map<locator *, locator_pvector, compare_locator> alignment_connection_map;

};

#endif
