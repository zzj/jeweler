#ifndef _EXON_NODE_H_
#define _EXON_NODE_H_

#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <Fasta.h>
#include "jeweler.hpp"
#include "transcript_info.hpp"
#include "transcript.hpp"
#include "laboratory/cigar_holder.hpp"
#include "landscape.plot.hpp"
using namespace BamTools;
using namespace std;

// The earrings always present as a pair, though sometimes you may
// only find one. 

// This class will investigate the possible allele specific pile up graph.



class ExonNode{
public :

	// start, end, and origin can determine the exon
	int start;
	int end;
	// 0 for non info
	// 1 for maternal
	// 2 for paternal
	int origin; 
	// only valid when origin != NO_INFO
	Exon * paried_exon;
	set<BamAlignment *> alignments;
	set<ExonNode *> in_nodes;
	set<ExonNode *> out_nodes;

};





#endif /* _EXON_NODE_H_ */
