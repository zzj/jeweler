
#ifndef _GRAPH_H_
#define _GRAPH_H_


#include "exon_node.hpp"
#include "path.hpp"
class Graph{
public:
	~Graph();
	
	set<ExonNode *> nodes;
	ExonNode * find_exon_node(int start, int end, int origin);
	ExonNode * add_exon_node(int start, int end, int origin, 
							 set<JewelerAlignment *> &reads,
							 set<JewelerAlignment *> &allele_reads);
	void add_edge(ExonNode * first, ExonNode *second, vector<JewelerAlignment *> reads);
	void add_edge(ExonNode * first, ExonNode *second, int num_reads);
	vector<ExonNode *> get_starting_nodes();
	void dump_graph(FILE * filename);
	int get_all_paths(vector<Path> &records);
	int traverse_graph( ExonNode * current_node, vector<ExonNode *> &path, vector<Path> &records, int origin);
	
};

#endif /* _GRAPH_H_ */
