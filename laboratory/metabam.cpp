#include "metabam.hpp"
#include "main/common.hpp"
#include "main/jeweler_alignment.hpp"

Metabam::Metabam() {
}

Metabam::Metabam(BamReader &reader) {
	initialize(reader);
}

void Metabam::initialize(BamReader &reader) {
	BamInfo::initialize(reader);
	total_num_alignments = 0;
	total_num_duplicate = 0;
	total_num_first_mate = 0;
	total_num_second_mate = 0;
	total_num_mapped = 0;
	total_num_mate_mapped = 0;
	total_num_paired = 0;
	total_num_proper_pair = 0;
	total_num_reverse_strand = 0;
	num_alignment_per_reference.resize(this->references.size(),0);
}

int Metabam::add_alignment(JewelerAlignment &al) {
	total_num_alignments++;
	num_alignment_per_reference[al.RefID]++;
	mapquality_list.push_back(al.MapQuality);
	if (al.IsDuplicate()) {
		total_num_duplicate++;
	}
	if (al.IsFirstMate()) {
		total_num_first_mate++;
	}
	if (al.IsSecondMate()) {
		total_num_second_mate++;
	}
	if (al.IsMapped()) {
		total_num_mapped++;
	}
	if (al.IsMateMapped()) {
		total_num_mate_mapped++;
	}
	if (al.IsProperPair()) {
		total_num_proper_pair++;
	}
	if (al.IsReverseStrand()) {
		total_num_reverse_strand++;
	}

	return 0;
}

void Metabam::dump_meta_data(FILE * log_file, bool is_readable) {
	size_t i;
	if (is_readable) {
		for (i = 0; i < references.size(); i++) {
			fprintf(log_file,"%s\t%d\n",references[i].RefName.c_str(),
					                    num_alignment_per_reference[i]);
		}
		fprintf(log_file,"%s\t%d\n","Total_alignments",total_num_alignments);
		fprintf(log_file,"%s\t%d\n","Duplicates",total_num_duplicate);
		fprintf(log_file,"%s\t%d\n","First_mate",total_num_first_mate);
		fprintf(log_file,"%s\t%d\n","Second_mate",total_num_second_mate);
		fprintf(log_file,"%s\t%d\n","Mapped",total_num_mapped);
		fprintf(log_file,"%s\t%d\n","Mate_mapped",total_num_mate_mapped);
		fprintf(log_file,"%s\t%d\n","Paired",total_num_paired);
		fprintf(log_file,"%s\t%d\n","Proper_pair",total_num_proper_pair);
		fprintf(log_file,"%s\t%d\n","Reverse_strand",total_num_reverse_strand);
	}
	else {
		;
	}
}

void Metabam::dump_mapquality_list(FILE *output) {
	for (size_t i=0; i<mapquality_list.size(); i++) {
		fprintf(output, "%d\n",mapquality_list[i]);
	}
}

