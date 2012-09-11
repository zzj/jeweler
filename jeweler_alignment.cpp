#include "jeweler_alignment.hpp"

void JewelerAlignment::dump_data(Jeweler::EarringsData::Read *read) {
    read->set_name(this->Name);
    for (auto i = genome_position.begin(); i != genome_position.end(); i++) {
        read->add_genome_position(i->second);
        read->add_read_position(i->first);
    }
    read->set_is_multiple_alignment(this->is_multiple_alignment);
    read->set_cigar_string(get_cigar_string(this->CigarData));
    return ;
}


void JewelerAlignment::set_is_multiple_alignment(bool a) {
    this->is_multiple_alignment = a;
}


int get_read_position(JewelerAlignment *al, int i) {
	if (al->IsReverseStrand()) {
		return (al->QueryBases.size() - i -1);
	}
	else{
		return i;
	}
}
