#ifndef _READ_MATCHER_CPP_
#define _READ_MATCHER_CPP_
#include <vector> 
#include "jeweler_alignment.hpp"

using namespace std;

class ReadMatcher {
public:
	vector<int> mismatch_read_locations;
	vector<int> mismatch_transcript_locations;
	vector<int> allele_read_locations;
	vector<int> allele_transcript_locations;
	vector<char> mismatchars;
	vector<char> mismatch_qualities;
	vector<int> transcript_aligned_locations;
	vector<int> read_aligned_locations;
	void output(JewelerAlignment *al);
};

#endif //  _READ_MISMATCHER_CPP_
