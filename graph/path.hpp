#ifndef _PATH_H_
#define _PATH_H_

#include "exon_node.hpp"
class Path{
public:
	Path(vector<ExonNode *> &);
	vector<ExonNode *> path;
	bool is_mirrored(Path & path);
	bool is_equal(Path & path);
	void dump_path(FILE * file);
	// contains allele specific exons or not
	bool is_informative();
	// contains NULL pointer or not
	bool is_valid();
// only contain one color nodes or not
bool is_compatible();

};
#endif
