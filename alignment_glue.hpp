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
#include "laboratory/cigar_holder.hpp"
using namespace BamTools;
using namespace std;


class AlignmentGlue{
public:
	map<int, vector<JewelerAlignment *> > refid2reads;
	map<string, vector<JewelerAlignment *> > name2reads;
	void glue_paired_alignments(JewelerAlignment *first, JewelerAlignment *second);
	int glue_cigar_data(JewelerAlignment *first, JewelerAlignment *second, int overlapped);	
	int get_skipped_region(const JewelerAlignment *al, int overlapped, 
						   vector<CigarOp> &cigar_data, int &skipped_length);
	

	int glue(vector<JewelerAlignment *> &in_reads, 
			 vector<JewelerAlignment *> &new_reads,
			 vector<JewelerAlignment *> &noused);
};

int output_bamalignment(JewelerAlignment *al);

#endif // _ALIGNMENT_GLUE_HPP_
