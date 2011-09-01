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
	vector<int> exon_start;
	vector<int> exon_end; //end is inclusive 
	vector<int> snp_pos;
	vector<char> alleles;
	transcript();
};

int load_gtf_file(string gtf_filename, vector<transcript> &transcripts);

#endif /* _GTF_H_ */
