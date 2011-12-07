#ifndef _LANDSCAPEPLOT_H_
#define _LANDSCAPEPLOT_H_

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

class LandscapePlot{
public:
	vector<int> unknown;
	vector<int> paternal;
	vector<int> maternal;
	vector<int> is_snp;
	vector<int> exon_jump;
	int num_maternal, num_paternal, num_unknown;
	string transcript_name;

	LandscapePlot(Transcript * maternal, 
				  Transcript * paternal, 
				  set<BamAlignment* > &unknown);
	
	int add_transcript_to_landscape(Transcript *,
									set<BamAlignment *> &reads,
									vector<int> & coverage);
								
	int generate_landscape_plot(FILE * info, FILE * output);

	int add_coverage(Transcript * transcript, BamAlignment *, vector<int>& coverage);
};



#endif /* _LANDSCAPE.PLOT_H_ */
