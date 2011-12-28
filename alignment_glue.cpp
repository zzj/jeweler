#include "alignment_glue.hpp"


int output_bamalignment(BamAlignment *al){
	fprintf(stdout,"%s\n",al->QueryBases.c_str());
	fprintf(stdout,"%d\n",al->Position);
	fprintf(stdout,"%s\n",get_cigar_string((*al)).c_str());
	return 0;
}

int AlignmentGlue::glue(vector<BamAlignment *> &in_reads, vector<BamAlignment *> &out_reads){
	int i;
	name2reads.clear();
	
	for ( i = 0; i < in_reads.size(); i++) {
		name2reads[ in_reads[i]->Name ].push_back(in_reads[ i ]);
	}

	for ( i = 0; i < in_reads.size(); i++) {
		if ( in_reads[i]->IsFirstMate() ) {
			if ( name2reads.find(in_reads[i]->Name) != name2reads.end()) {
				if ( name2reads[in_reads[i]->Name].size() > 2 ) {
					fprintf(stdout, 
							"Oops, there are mulitple alignments within a gene?!\n");
					vector<BamAlignment *> alignments=
						name2reads[in_reads[i]->Name];
					for (int i = 0; i < alignments.size(); i++){
						output_bamalignment(alignments[i]);
					}
				}
			}
		}
	}
	
	return 0;
}
