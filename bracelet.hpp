#ifndef _BRACELET_HPP_
#define _BRACELET_HPP_
#include "transcript_info.hpp"
#include <vector> 
#include <cstdio>
#include <algorithm>
using namespace std;
class Bracelet{
public:
	vector<vector<string> > reads;
	vector<vector<int> > results;
	vector<vector<int> > related_transcripts;
	vector<TranscriptInfo *> info;
	
	Bracelet(vector<TranscriptInfo *> info);
	int intersect(vector<string> &a, vector<string> &b);
	int analyze();
	int dump(FILE *file);

};
#endif
