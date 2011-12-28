#ifndef _ALIGNMENT_GLUE_HPP_
#define _ALIGNMENT_GLUE_HPP_
#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <map>
#include <Fasta.h>
#include "laboratory/cigar_holder.hpp"
using namespace BamTools;
using namespace std;


class AlignmentGlue{
public:
	map<int, vector<BamAlignment *> > refid2reads;
	map<string, vector<BamAlignment *> > name2reads;
	int glue(vector<BamAlignment *> &in_reads, vector<BamAlignment *> &out_reads);
};

#endif // _ALIGNMENT_GLUE_HPP_
