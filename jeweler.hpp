#ifndef _JEWELER_H_
#define _JEWELER_H_

#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <map>
#include <set>
#include "fasta.hpp"
#include "common.hpp"
#include "gtf.hpp"
#include "rna_read.hpp"
#include "transcript_info.hpp"
#include "transcript.hpp"
#include "earrings.hpp"
#include "laboratory/sewing_machine.hpp"

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <Fasta.h>

using namespace BamTools;

class jeweler{
public:
	string info_filename;
	string mamf_filename;
	FILE * log_file;
	vector<TranscriptInfo *> transcripts_info;

	SewingMachine *sm;
	jeweler(int argc, char *argv[]);
	// function to run the analysis 
	int run();

	// load info file, which is generated by split_gtf.php
	// each line in info file is a cufflinks gene.
	int load_info_file();
	int load_mamf_file();

};

#endif /* _JEWELER_H_ */
