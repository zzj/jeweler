
#ifndef _RNA_READ_H_
#define _RNA_READ_H_
#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <vector>
#include <boost/algorithm/string.hpp>
#include "common.hpp"
using namespace std;

class rna_read_key{
public:
	rna_read_key(string tname, string rname);
	string transcript_name;
	string rna_read_name;
};

bool operator < (const rna_read_key &a, const rna_read_key &b);

class rna_read_query{
public:
	string seq;
	string name;
	string target; // transcript id;
	string flag_field;
	int size;
	vector<int> target_start;
	vector<int> query_start;
	vector<int> block_size;

	int first_end;       // the rightmost aligned position of the left end read
	int second_start;    // the leftmost aligned position of the right end read
	
	int start_in_query;
	int end_in_query;
	int mismatch;
	int matches;
	int target_gap_num;
	int target_gap_size;
	int rep_match;
	int unknown_size;
	int query_gap_num;
	int query_gap_size;
	// the source of the reads, 
	// 0 means unknown
	// 1 means paternal
	// 2 means maternal
	int source_id; 
	bool is_ignored;
	bool is_initialized;
	bool is_merged;

	rna_read_query();
	int set_seq(string seq);
	bool is_reversed();
};

bool operator < (const rna_read_query& a, const rna_read_query&  b);
bool is_better_alignment(const rna_read_query* a, const rna_read_query*  b);
void recover_original_read(string &seq);
int load_psl_file(string psl_filename, vector<rna_read_query *> & queries);

#endif /* _RNA_READ_H_ */
