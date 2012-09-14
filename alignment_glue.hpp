#ifndef _ALIGNMENT_GLUE_HPP_
#define _ALIGNMENT_GLUE_HPP_

#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <map>
#include <set>
#include <Fasta.h>
#include <algorithm>
using namespace BamTools;
using namespace std;

class JewelerAlignment;

class AlignmentGlue{
public:
	void glue_paired_alignments(JewelerAlignment *first, JewelerAlignment *second);
	int glue_cigar_data(JewelerAlignment *first, JewelerAlignment *second, int overlapped);	

	int glue(vector<JewelerAlignment *> &in_reads, 
			 vector<JewelerAlignment *> &new_reads,
			 vector<JewelerAlignment *> &noused);

private:
	map<int, vector<JewelerAlignment *> > refid2reads;
	map<string, vector<JewelerAlignment *> > name2reads;

};

int output_bamalignment(JewelerAlignment *al);

#endif // _ALIGNMENT_GLUE_HPP_
