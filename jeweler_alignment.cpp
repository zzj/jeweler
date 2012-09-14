#include "jeweler_alignment.hpp"
#include "laboratory/cigar_holder.hpp"
#include "laboratory/sewing_machine.hpp"

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


int get_read_position(JewelerAlignment *al, const int i) {
	if (al->IsReverseStrand()) {
		return (al->QueryBases.size() - i -1);
	}
	else{
		return i;
	}
}

void JewelerAlignment::set_is_multiple_alignment(const SewingMachine *sm) {
    if (sm->is_multiple_alignment(this->Name)) {
        this->set_is_multiple_alignment(true);
    }
    else {
        this->set_is_multiple_alignment(false);
    }
    return ;
}

void JewelerAlignment::jeweler_initialize(const SewingMachine *sm) {
    assert(sm);
    

}
