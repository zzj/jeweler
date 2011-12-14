#ifndef _EXON_NODE_H_
#define _EXON_NODE_H_

#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <Fasta.h>
#include <set> 
#include "../constants.hpp"
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

	ExonNode * paried_exon;
	set<BamAlignment *> alignments;

	vector<ExonNode *> in_nodes;
	set<ExonNode *> in_nodes_checklist;
	vector<int> in_edge_weight;
	vector<set<BamAlignment *> > in_edge_reads;
	
	vector<ExonNode *> out_nodes;
	set<ExonNode *> out_nodes_checklist;
	vector<int> out_edge_weight;
	vector<set<BamAlignment *> > out_edge_reads;

	ExonNode(int start, int end, int origin);
	string detach();
	bool is_mirrored(ExonNode *a);
	bool is_equal(ExonNode *a);
};





#endif /* _EXON_NODE_H_ */
