
#include "sewing_machine.hpp"

locator::locator(const JewelerAlignment &al) {
	this->ref_id = al.RefID;
	this->position = al.Position;
	this->cigar_string = get_cigar_string(al);
	al.GetTag("NM", this->edit_distance);
}

bool compare_locator::operator()(const locator *a, const locator *b) {
	if (a->ref_id!=b->ref_id) return a->ref_id<b->ref_id;
	else return a->position<b->position;
}

SewingMachine::SewingMachine() {
}

SewingMachine::SewingMachine(BamReader &reader) {
	this->initialize(reader);
}

void SewingMachine::initialize(BamReader &reader) {
	BamInfo::initialize(reader);
	this->num_multiple_reads_per_chr.resize(references.size(),0);
}

string SewingMachine::get_reference_sequence(FastaReference &fr, JewelerAlignment &al) {
	string ret;
	// TODO : read from reference genome to get the mismatch position
	// or other useful information.
	return ret;
}

string SewingMachine::get_full_name(const JewelerAlignment &al) {
	string is_first;
	if (al.IsFirstMate()) {
		is_first = "\t1";
	}
	else {
		is_first = "\t2";
	}
    return al.Name + is_first;
}

int SewingMachine::count_alignment(const JewelerAlignment &al) {
    string full_name = get_full_name(al);
    if ( multiple_alignment_map.find(full_name) ==
         multiple_alignment_map.end()) {
        multiple_alignment_map[full_name] = 1;
    }
    else {
        multiple_alignment_map[full_name] ++;
    }
	return 0;
}

int SewingMachine::add_alignment(const JewelerAlignment &al) {
	locator *l =new locator(al);
    string full_name = get_full_name(al);
    if (multiple_alignment_map[full_name] > 1)
        seqs[full_name].push_back(l);
	return 0;
}

int SewingMachine::output_alignment_map_core(FILE * file, bool only_mulitple_reads) {
	auto i=seqs.begin();
	size_t j;
	for (; i!=seqs.end();i++) {
		if (only_mulitple_reads && i->second.size()<2) continue;
		fprintf(file,"%s\t%zd\t",i->first.c_str(),i->second.size());
		for (j=0;j<i->second.size();j++) {
			fprintf(file,"%s\t%d\t%s\t%d\t",
					references[i->second[j]->ref_id].RefName.c_str(),
					i->second[j]->position,
					i->second[j]->cigar_string.c_str(),
					i->second[j]->edit_distance);
		}
		fprintf(file,"\n");
	}
	return 0;
}

void SewingMachine::load_multiple_alignments_set(FILE *file ) {
	char refname[1000], cigar[1000];
	char read_id[1000];
	char is_first[100];
	int pos, distance, size;

	int i;
	int last;
	while((last = fscanf(file, "%s%s%d", read_id, is_first, &size)) == 3) {
		// ignored
		for (i = 0; i < size; i++) {
			fscanf(file, "%s%d%s%d", refname, &pos, cigar, &distance);
		}

		// only record the multiple alignment
		if (size >= 2) {
			multiple_alignment_set.insert(read_id);
		}
	}
}

int SewingMachine::output_alignment_map(FILE * file) {
	return output_alignment_map_core(file, /*only output multiple reads?*/ false);
}

int SewingMachine::output_multiple_alignment_map(FILE * file) {
	return output_alignment_map_core(file, /*only output multiple reads?*/ true);
}

int SewingMachine::build_alignment_connection_map_core() {
	auto i = seqs.begin();
	size_t j,k;
	alignment_connection_map.clear();
	for ( ; i != seqs.end(); i++) {
		if (i->second.size() < 2) continue;
		for (j = 0; j < i->second.size(); j++) {
			// statistics
			num_multiple_reads_per_chr[i->second[j]->ref_id]++;

			// map
			for (k = 0;k < i->second.size(); k++) {
				if (j != k)
					alignment_connection_map[i->second[j]].push_back(i->second[k]);
			}
		}
	}
	return 0;
}


int SewingMachine::output_alignment_connection_map_core(FILE * file) {
	auto i=alignment_connection_map.begin();
	for (; i!=alignment_connection_map.end(); i++) {
		size_t j=0;
		for (; j!=i->second.size(); j++) {
			fprintf(file, "%s\t%d\t%s\t%d\n",
					references[i->first->ref_id].RefName.c_str(),
					i->first->position,
					references[i->second[j]->ref_id].RefName.c_str(),
					i->second[j]->position);
		}
	}
	return 0;
}

int SewingMachine::output_alignment_connection_map(FILE * file) {
	build_alignment_connection_map_core();
	output_alignment_connection_map_core(file);
	return 0;
}
