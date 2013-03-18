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
	PileupPlot(Transcript * maternal,
			   Transcript * paternal,
			   set<JewelerAlignment *> &unknown,
			   set<JewelerAlignment *> &multiple);

	void generate_pileup_plot(FILE * info, FILE * output);

private:
	int add_transcript_to_pileup_filter(Transcript * transcript,
										set<JewelerAlignment *> &reads,
										vector<int> & coverage,
										set<JewelerAlignment *> &multiple);

	int add_transcript_to_pileup(Transcript *,
                                 set<JewelerAlignment *> &reads,
                                 vector<int> & coverage);

	int add_coverage(Transcript * transcript, JewelerAlignment *, vector<int>& coverage);

private:
	vector<int> unknown;
	vector<int> paternal;
	vector<int> maternal;
	vector<int> multiple;
	vector<int> is_snp;
	vector<int> exon_jump;
	vector<int> genome_pos;
	int num_maternal, num_paternal, num_unknown, num_exons;
	string transcript_id;
};

#endif /* _PILEUP.PLOT_H_ */
