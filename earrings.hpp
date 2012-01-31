#ifndef _EARRINGS_H_
#define _EARRINGS_H_



#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <Fasta.h>
#include "jeweler.hpp"
#include "transcript_info.hpp"
#include "transcript.hpp"
#include "transcript_mismatcher.hpp"
#include "laboratory/cigar_holder.hpp"
#include "laboratory/sewing_machine.hpp"
#include "pileup.plot.hpp"
#include "alignment_glue.hpp"
using namespace BamTools;
using namespace std;

// The earrings always present as a pair, though sometimes you may
// only find one. 

// This class will investigate the possible allele specific pile up graph.


class Earrings{
public :
	Earrings(TranscriptInfo *, SewingMachine * sm);

	~Earrings();
	TranscriptInfo * info;
	TranscriptMismatcher *mismatcher;
	SewingMachine *sm;
	vector<Transcript *> maternal_transcripts;
	vector<Transcript *> paternal_transcripts;
	vector<JewelerAlignment *> bam_reads;
	vector<JewelerAlignment *> unaligned_reads;
	vector<JewelerAlignment *> noused_reads;
	vector<JewelerAlignment *> compatible_reads;
	

	int num_total_reads;
	// find the reads that are mulitple aligned to other places
	int count_multiple_alignments();
	// create paternal and maternal transcripts databases
	// annotate SNPs for paternal  and maternal transcripts' sequence
	int load_transcript_data(TranscriptInfo * ti); // load transcript data 

	// load cufflinks' gtf file 
	int transcript_helper(string read_file,string gtf_file, 
						  vector<Transcript *> &transcripts);

	// load JewelerAlignment into bam_reads
	int load_read_data(TranscriptInfo *ti);
	// align reads to maternal or paternal transcripts
	int align_reads();

	// test allele specific transcript
	int test_allele_specific_transcript();

	// build graph
	int build_graph();
	
	// study compatible reads
	int study_compatible_reads();

	// get compatible reads
	int get_compatible_reads(vector<set<JewelerAlignment*> >& read_lists );

	void test_memory_leak();
};


#endif /* _EARRINGS_H_ */
