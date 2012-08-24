#ifndef _BRACELET_HPP_
#define _BRACELET_HPP_
#include "jeweler.hpp"
#include <vector> 
#include <cstdio>
#include <algorithm>
#include <set>
using namespace std;
class JewelerInfo;
class Bracelet{
public:
	vector<vector<string> > reads;
	vector<set<string> > reads_index;
	vector<vector<int> > results;
	vector<vector<int> > related_transcripts;
	JewelerInfo *jeweler_info;
	string gene_id;
	string result_folder;
	Bracelet(JewelerInfo *jeweler_info);
	int intersect(vector<string> &a, vector<string> &b);
	int analyze();
	int dump(FILE *file);
    void dump_shared_pileup(FILE * fd,
                           int original_id,
                           int target_id);

	int dump(FILE * file, string root);

};
#endif
