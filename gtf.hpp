#ifndef _GTF_H_
#define _GTF_H_
#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <vector>
#include <boost/algorithm/string.hpp>
#include "common.hpp"

using namespace std;
class transcript{
public:
	string name;
	string seq;
	string chr;  // storing chromosome name
	vector<int> exon_start;
	vector<int> exon_end; //end is inclusive 
	vector<int> snp_pos;
	vector<char> alleles;
	vector<int> noninformative_mismatches; // number of non informative mismatches per base
	vector<int> Anoninformative_mismatches; // number of non informative mismatches equaling to A per base
	vector<int> Cnoninformative_mismatches; // number of non informative mismatches equaling to C per base
	vector<int> Gnoninformative_mismatches; // number of non informative mismatches equaling to G per base
	vector<int> Tnoninformative_mismatches; // number of non informative mismatches equaling to T per base
	vector<int> genome_pos;  // storing the within chromosome pos of each base in the transcript
	transcript();
};

int load_gtf_file(string gtf_filename, vector<transcript *> &transcripts);

#endif /* _GTF_H_ */
