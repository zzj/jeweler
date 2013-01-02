#include "alignment_glue.hpp"
#include "jeweler_alignment.hpp"
#include "laboratory/cigar_holder.hpp"

int output_bamalignment(JewelerAlignment *al) {
	fprintf(stdout,"%s\n",al->Name.c_str());
	if (al->IsFirstMate()) {
		fprintf(stdout,"first one\n");
	}
	else {
		fprintf(stdout,"second one\n");
	}
	fprintf(stdout,"%s\n",al->QueryBases.c_str());
	fprintf(stdout,"%s\n",al->Qualities.c_str());
	fprintf(stdout,"%d\n",al->Position + 1);
	fprintf(stdout,"%s\n",get_cigar_string((*al)).c_str());
	fprintf(stdout, "\n");
	return 0;
}


void AlignmentGlue::glue_paired_alignments(JewelerAlignment *first, JewelerAlignment *second) {
    first->glue(second);
}


int AlignmentGlue::glue(vector<JewelerAlignment *> &in_reads,
                        vector<JewelerAlignment *> &new_reads,
                        vector<JewelerAlignment *> &noused_reads) {
    size_t i;

    int ignore_unpaired_read = 0;
    int ignore_unproper_read = 0;
    int num_first_mate = 0;

    set<JewelerAlignment *> checklist;

    name2reads.clear();

    for (i = 0; i < in_reads.size(); i++) {
        name2reads[ in_reads[i]->Name ].push_back(in_reads[ i ]);
    }

    for (i = 0; i < in_reads.size(); i++) {
        if (in_reads[i]->IsFirstMate()) {
            num_first_mate ++;
            if (name2reads.find(in_reads[i]->Name) != name2reads.end()) {
                vector<JewelerAlignment *> alignments = name2reads[in_reads[i]->Name];
                if (alignments.size() == 1) {
                    // TODO:
                    // not paired, but still counted as a valid
                    // alignment.
                    ignore_unpaired_read ++;
                    continue;
                }
                else if (alignments.size()==2) {
                    if ((alignments[1]->IsFirstMate() &&
                         alignments[0]->IsSecondMate()) ||
                        (alignments[0]->IsFirstMate() &&
                         alignments[1]->IsSecondMate())) {
                        if (alignments[0]->Position < alignments[1]->Position) {
                            glue_paired_alignments(alignments[0], alignments[1]);
                            new_reads.push_back(alignments[0]);
                            checklist.insert(alignments[0]);
                        }
                        else if (alignments[1]->Position >= alignments[0]->Position) {
                    		glue_paired_alignments(alignments[1], alignments[0]);
							new_reads.push_back(alignments[1]);
							checklist.insert(alignments[1]);
						}
					}
					else {
						//fprintf(stdout, "not properly mapped\n");
						ignore_unproper_read ++;
					}
				}
				else {
					//fprintf(stdout,
					//"Oops, there are mulitple alignments within a gene?!\n");
					//for (int i = 0; i < alignments.size(); i++) {
					//output_bamalignment(alignments[i]);
					//}
				}
			}
		}
	}
	//fprintf(stdout, " %d unpaired and %d not propered and %d of first paried\n", ignore_unpaired_read, ignore_unproper_read, num_first_mate);
	for (i = 0; i < in_reads.size(); i++) {
		if (checklist.find(in_reads[i]) == checklist.end()) {
			noused_reads.push_back(in_reads[i]);
			continue;
			// push into single paired reads

			// new_reads.push_back(in_reads[i]);
			// for (j = 0; j < in_reads[i]->QueryBases.size(); j ++) {
			// 	in_reads[i]->read_position.push_back(get_read_position(in_reads[i],j));
			// }

		}
	}

	return 0;
}
