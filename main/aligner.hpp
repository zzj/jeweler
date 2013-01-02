#ifndef _ALIGNER_HPP_
#define _ALIGNER_HPP_
#include <cstring>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>  
using namespace std;

class TranscriptomeAligner{
public:
	string left_unmapped_file, right_unmapped_file;
	vector<string> ref_file;
	TranscriptomeAligner();
	void run_command(string command);
	void align(string &left_unmapped_file, string &right_unmapped_file, 
			   vector<string>  &reference_sequences, vector<string> &ids,
			   string &merged_fasta_file,
			   string &output_file);
};

#endif
