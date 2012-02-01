#ifndef _TRANSCRIPT_MISMATCHER_H_
#define _TRANSCRIPT_MISMATCHER_H_

#include "transcript.hpp"
#include "constants.hpp"
#include <boost/dynamic_bitset.hpp>
#include <boost/math/distributions/poisson.hpp>
#include "transcript_info.hpp"
#include <algorithm> 
#include <functional>
#include <numeric>
#include <cmath>

using namespace std;

class TranscriptMismatcher{
	
public:

	TranscriptMismatcher();
	
	int initialize();
	
	int add_transcript(Transcript *, int origin);
	
	int add_mismatches(Transcript *transcript, JewelerAlignment *al, 
					   vector<int> &transcript_locations,
					   vector<int> &read_locations,
					   vector<int> &transcript_mismatch_locations,
					   vector<char> &read_mismatch_qualities,
					   vector<char> & mismatchars);

	int dump(FILE *);
	int write(FILE *);
	
	map<int, int> genome_pos2idx;
	map<int, char> genome_pos2paternal;
	map<int, char> genome_pos2maternal;
	
	vector<int> mismatches; // number of non informative mismatches
							// per base 
	set<JewelerAlignment *> reads;
	vector<int> coverage;
	// the list of mismatching reads per base pair
	vector<set<JewelerAlignment *> > mismatched_reads;
	// the list of mismatching positions 
	vector<vector<char > > read_mismatch_qualities;
	vector<vector<char > > read_qualities;

	// number of non informative mismatches equaling to A per base
	vector<int> A_mismatches; 
	// number of non informative mismatches equaling to C per base
	vector<int> C_mismatches; 
	// number of non informative mismatches equaling to G per base
	vector<int> G_mismatches; 
	// number of non informative mismatches equaling to T per base
	vector<int> T_mismatches; 
	// number of non informative mismatches equaling to N per base
	vector<int> N_mismatches; 
	
};

class TranscriptMismatcherAnalyzer{
public:
	vector<int> genome_locations;
	vector<int> coverages;
	vector<int> mismatches;
	vector<char> maternal_seq, paternal_seq;
	vector<vector<char> > read_mismatch_qualities;
	vector<vector<char> > read_qualities;
	
	// number of non informative mismatches equaling to A per base
	vector<int> A_mismatches; 
	// number of non informative mismatches equaling to C per base
	vector<int> C_mismatches; 
	// number of non informative mismatches equaling to G per base
	vector<int> G_mismatches; 
	// number of non informative mismatches equaling to T per base
	vector<int> T_mismatches; 
	// number of non informative mismatches equaling to N per base
	vector<int> N_mismatches; 
	int num_reads;
	int num_locations;
	boost::dynamic_bitset<> is_consistent_mismatches;
	map<char, double> error_rate;
	vector<double> p_values;

	bool is_initialized;
	
	TranscriptMismatcherAnalyzer();
	TranscriptMismatcherAnalyzer(vector<TranscriptInfo *> &ti);

	int append(FILE *);
	
	int end_loading();
	
	int calculate_error();

	int calculate_p_value();

	int mark_consistent_mismatch();

	int analyze();
};


#endif
