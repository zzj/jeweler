#ifndef _PILEUPPLOT_H_
#define _PILEUPPLOT_H_

#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <vector>
#include <algorithm>
#include <set>
#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;

#include "common.hpp"
#include "transcript.hpp"
using namespace std;

class PileupPlot{
public:
	vector<int> unknown;
	vector<int> paternal;
	vector<int> maternal;
	vector<int> is_snp;
	vector<int> exon_jump;
	vector<int> genome_pos;
	int num_maternal, num_paternal, num_unknown;
	string transcript_name;

	PileupPlot(Transcript * maternal, 
				  Transcript * paternal, 
				  set<JewelerAlignment* > &unknown);
	
	int add_transcript_to_pileup(Transcript *,
									set<JewelerAlignment *> &reads,
									vector<int> & coverage);
								
	int generate_pileup_plot(FILE * info, FILE * output);

	int add_coverage(Transcript * transcript, JewelerAlignment *, vector<int>& coverage);
};



#endif /* _PILEUP.PLOT_H_ */
