
#ifndef _GRAPH_H_
#define _GRAPH_H_


#include "exon_node.hpp"
#include "path.hpp"
class Graph{
public:
	set<ExonNode *> nodes;
	ExonNode * find_exon_node(int start, int end, int origin);
	ExonNode * add_exon_node(int start, int end, int origin);
	int add_edge(ExonNode * first, ExonNode *second, vector<BamAlignment *> reads);
	int add_edge(ExonNode * first, ExonNode *second, int num_reads);
	vector<ExonNode *> get_starting_nodes();
	int dump_graph(FILE * filename);
	int get_all_paths(vector<Path> &records);
	int traverse_graph( ExonNode * current_node, vector<ExonNode *> &path, vector<Path> &records);
};

#endif /* _GRAPH_H_ */
