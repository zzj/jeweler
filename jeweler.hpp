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
#include "gtf.hpp"
#include "rna_read.hpp"
#include "transcript.hpp"
#include "earrings.hpp"
#include "bracelet.hpp"
#include "transcript_mismatcher.hpp"
#include "laboratory/sewing_machine.hpp"

#include <api/BamReader.h>
#include <api/BamWriter.h>

#include <Fasta.h>

using namespace BamTools;


#include <boost/filesystem.hpp>

using boost::filesystem::path;
using boost::filesystem::create_directory;


class JewelerInfo{
public:
	string maternal_strain_ref_file;
	string paternal_strain_ref_file;
	string maternal_strain_id;
	string paternal_strain_id;
	string bam_file;
	string alias;
	string result_file;
	string result_folder;
	string gtf_input_file;

	BamReader bam_reader;
	RefVector references;
	FastaReference *maternal_fasta, *paternal_fasta;

	vector<Transcript *> transcripts;

	map<string, vector<Transcript *> > gene_id2transcripts;
	vector<string> gene_id;
	JewelerInfo(int argc, char * argv[]);
	~JewelerInfo();
	int check_args(int &i, char *  argv[], const char *name, string &a);
	int build_gene_id2transcripts();
	int get_refID(string chr);
};

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

	bool is_earrings;
	bool is_bracelet;
	bool is_mismatch_analyzer;
	// function to run the analysis 
	int run();

	int load_mamf_file();
	
};

#endif /* _JEWELER_H_ */
