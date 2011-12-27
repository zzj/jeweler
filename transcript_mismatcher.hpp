#ifndef _TRANSCRIPT_MISMATCHER_H_
#define _TRANSCRIPT_MISMATCHER_H_

#include "transcript.hpp"
#include "constants.hpp"
class TranscriptMismatcher{
	
public:

	TranscriptMismatcher();
	
	int initialize();
	
	int add_transcript(Transcript *, int origin);
	
	int add_mismatches(Transcript *transcript, BamAlignment *al, 
					   vector<int> &locations,
					   vector<int> &transcript_locations,
					   vector<char> & mismatchars);

	int dump(FILE *);
	
	map<int, int> genome_pos2idx;
	map<int, char> genome_pos2paternal;
	map<int, char> genome_pos2maternal;
	
	vector<int> mismatches; // number of non informative mismatches
							// per base
	set<BamAlignment *> reads;
	vector<int> coverage;
	// the list of mismatching reads per base pair
	vector<set<BamAlignment *> > mismatched_reads;

	// number of non informative mismatches equaling to A per base
	vector<int> A_mismatches; 
	// number of non informative mismatches equaling to C per base
	vector<int> C_mismatches; 
	// number of non informative mismatches equaling to G per base
	vector<int> G_mismatches; 
	// number of non informative mismatches equaling to T per base
	vector<int> T_mismatches; 
	
};

#endif
