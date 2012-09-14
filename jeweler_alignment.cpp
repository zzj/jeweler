#include "jeweler_alignment.hpp"
#include "laboratory/cigar_holder.hpp"
#include "laboratory/sewing_machine.hpp"


void AlignmentExpert::study_matched_seq(JewelerAlignment *al,
                                        int genome_start,
                                        int alignment_start,
                                        int length) {
    return ;
}

void AlignmentExpert::study_only_read_seq(JewelerAlignment *al,
                                          int genome_start,
                                          int alignment_start,
                                          int length) {
    return ;
}

void AlignmentExpert::study_only_genome_seq(JewelerAlignment *al,
                                            int genome_start,
                                            int alignment_start,
                                            int length) {
    return ;
}


void AlignmentExpert::study_neither_exist_seq(JewelerAlignment *al,
                                               int genome_start,
                                               int alignment_start,
                                               int length) {
    return ;
}


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
    this->set_is_multiple_alignment(sm);
}


int JewelerAlignment::GenomeStartPosition() {
    return this->Position + 1;
}

void JewelerAlignment::investigate(AlignmentExpert *ae) {
    // must perform the compatible test first!
    // justify whether the sequences contains the JewelerAlignment

    // not sure this is the case, but the coordinates are messed up
    // sometimes. TODO: read the samtools's specification, and make
    // the coordinates right.
    int genome_start=this->GenomeStartPosition();
    int alignment_start=0;

    const std::vector< CigarOp > &cigar_data = this->CigarData;
    vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
    vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();

    for (; cigar_iter != cigar_end; ++cigar_iter) {
        const CigarOp& op = (*cigar_iter);

        switch (op.Type) {

            // for 'M', '=', 'X' - aligned string
        case (Constants::BAM_CIGAR_MATCH_CHAR)    :
        case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
        case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
            // the beginning and end of the matched sequence must
            // be belong to the same sequences. <
            ae->study_matched_seq(this, genome_start, alignment_start, op.Length);
            genome_start += op.Length;
            alignment_start += op.Length;
            break;
            // none aligned string in reads
        case (Constants::BAM_CIGAR_INS_CHAR)      :
        case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
            ae->study_only_read_seq(this, genome_start, alignment_start, op.Length);
            alignment_start += op.Length;
            break;

            // none aligned string in reference genome
        case 'J'      :
        case (Constants::BAM_CIGAR_DEL_CHAR) :
        case (Constants::BAM_CIGAR_PAD_CHAR) :
        case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
            ae->study_only_genome_seq(this, genome_start, alignment_start, op.Length);
            genome_start += op.Length;
            break;

            // for 'H' - hard clip, do nothing to AlignedBases, move to next op
        case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
            ae->study_neither_exist_seq(this, genome_start, alignment_start, op.Length);
            break;

            // invalid CIGAR op-code
        default:
            const string message = string("invalid CIGAR operation type: ") + op.Type;
            fprintf(stdout, "%s\n", message.c_str());
            exit(0);
        }
    }
	return ;

}
