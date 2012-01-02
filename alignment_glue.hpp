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
	map<int, vector<BamAlignment *> > refid2reads;
	map<string, vector<BamAlignment *> > name2reads;
	int glue_paired_alignments(BamAlignment *first, BamAlignment *second);
	int glue_cigar_data(BamAlignment *first, BamAlignment *second, int overlapped);	
	int get_skipped_region(BamAlignment *al, int overlapped, 
						   vector<CigarOp> &cigar_data, int &skipped_length);
	

	int glue(vector<BamAlignment *> &in_reads, 
			 vector<BamAlignment *> &new_reads,
			 vector<BamAlignment *> &noused);
};

int output_bamalignment(BamAlignment *al);

int cigar_trim(BamAlignment *al);
#endif // _ALIGNMENT_GLUE_HPP_
