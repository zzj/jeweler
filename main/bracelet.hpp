#ifndef _BRACELET_HPP_
#define _BRACELET_HPP_
#include <vector> 
#include <cstdio>
#include <algorithm>
#include <set>
#include <fstream>
#include <memory>
#include "proto/jeweler.pb.h"

using namespace std;
class JewelerInfo;
class ZMegaFile;

class Bracelet{
public:
	Bracelet(JewelerInfo *jeweler_info, int num_test_case = -1);
	int intersect(vector<string> &a, vector<string> &b);
	int analyze();
    void dump_shared_pileup(Jeweler::BraceletData::RelatedTranscript *fd,
                           int original_id,
                           int target_id);

	int dump(fstream *, string root);
    int dump_all();
private:
    Jeweler::BraceletData *data;
	vector<vector<string> > reads;
	vector<map<string, int> > reads_index;
	vector<vector<int> > results;
	vector<vector<int> > related_transcripts;
	JewelerInfo *jeweler_info;
	string gene_id;
	string result_folder;
    shared_ptr<ZMegaFile> zmf;
};

double get_coverage_rate(const Jeweler::EarringsData::Read &origin,
                         const map<int, int> &coverage);

void add_coverage_details(const Jeweler::EarringsData::Read &origin,
                          const Jeweler::EarringsData::Read &target,
                          map<int, map<int, int> > &genome_position_map);

void add_coverage(const Jeweler::EarringsData::Read &origin,
                          map<int, int> &coverage);
#endif
