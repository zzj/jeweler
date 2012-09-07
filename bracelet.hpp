#ifndef _BRACELET_HPP_
#define _BRACELET_HPP_
#include <vector> 
#include <cstdio>
#include <algorithm>
#include <set>
#include <fstream>
#include "proto/jeweler.pb.h"

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
    void dump_shared_pileup(Jeweler::BraceletData::RelatedTranscript *fd,
                           int original_id,
                           int target_id);

	int dump(fstream *, string root);
    int dump_all();
private:
    Jeweler::BraceletData *data;
};
#endif
