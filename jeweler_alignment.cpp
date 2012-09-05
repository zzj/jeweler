#include "jeweler_alignment.hpp"

void JewelerAlignment::dump_data(Jeweler::EarringsData::Read *read) {
    read->set_name(this->Name);
    read->set_genome_position(this->Position);
    read->set_cigar_string(get_cigar_string(this->CigarData));
    return ;
}

int get_read_position(JewelerAlignment *al, int i) {
	if (al->IsReverseStrand()) {
		return (al->QueryBases.size() - i -1);
	}
	else{
		return i;
	}
}
