#include "jeweler_alignment.hpp"
#include "laboratory/cigar_holder.hpp"
#include "laboratory/sewing_machine.hpp"
#include "constants.hpp"

void AlignmentExpert::initialize(JewelerAlignment *al) {
    return ;
}

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

void AlignmentExpertStarter::initialize(JewelerAlignment *al) {
    al->genome_position.clear();
}

void AlignmentExpertStarter::study_matched_seq(JewelerAlignment *al,
                                               int genome_start,
                                               int alignment_start,
                                               int length) {
    for (int i = 0; i < length; i ++ ){
        al->genome_position.push_back(genome_start + i);
    }
    return ;
}

void AlignmentExpertStarter::study_only_read_seq(JewelerAlignment *al,
                                          int genome_start,
                                          int alignment_start,
                                          int length) {
    for (int i = 0; i < length; i ++ ){
        al->genome_position.push_back(NOT_FOUND);
    }
    return ;
}

void JewelerAlignment::dump_data(Jeweler::EarringsData::Read *read) {
    assert(this->QueryBases.size() == this->genome_position.size());
    read->set_name(this->Name);
    for (auto i = genome_position.begin(); i != genome_position.end(); i++) {
        read->add_genome_position(*i);
    }
    read->set_seq(this->QueryBases);
    read->set_is_multiple_alignment(this->is_multiple_alignment);
    read->set_cigar_string(get_cigar_string(this->CigarData));
    read->set_is_second_truncated(this->is_second_truncated);
    read->set_glue_position(this->glue_position);
    read->set_tail_length(this->tail_length);
    read->set_head_length(this->head_length);
    read->set_is_reverse_strand(this->IsReverseStrand());
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
    AlignmentExpertStarter aes;
    this->investigate(&aes);
}


int JewelerAlignment::get_start_position() {
    return this->Position + 1;
}

int JewelerAlignment::get_end_position() {
    return this->GetEndPosition() + 1;
}


void JewelerAlignment::glue(JewelerAlignment *that) {
	// DEBUG:
	// output_bamalignment(first);
	// output_bamalignment(second);
	size_t i, j;
	// int first_start	 = this->Position + 1;
	size_t first_end	 = this->get_end_position();
	size_t second_start = that->get_start_position();
	size_t second_end	 = that->get_end_position() ;
	int overlapped	 = second_start - first_end;
	// DEBUG:
	// Check weather the first alignment is always before the second
	// alignment
    this->is_second_truncated = this->IsFirstMate();
    this->glue_position = this->Length;
    this->tail_length = that->Length;
    this->head_length = this->Length;

	if (overlapped>=0) {
		// Great news! no overlapped region between first and second
		// reads.
		this->Length += that->Length;
		this->Qualities += that->Qualities;
		string pad(overlapped, 'N');
		this->QueryBases += that->QueryBases;
		if (overlapped > 0) {
			this->CigarData.push_back(CigarOp('J',overlapped));
		}

		this->CigarData.insert(this->CigarData.end(),
								that->CigarData.begin(), that->CigarData.end());
	}
	else if (second_end > first_end) {
		// bad news! Overlapped region need to be removed from the
		// second read.
		int skipped_read_length = 0 , skipped_genome_length = -overlapped;
        vector<CigarOp> new_cigar_data;
        that->get_skipped_region(skipped_genome_length,
                                 new_cigar_data, skipped_read_length);
        this->Length += that->Length - skipped_read_length;
        this->QueryBases += that->QueryBases.substr(skipped_read_length);
        this->Qualities += that->Qualities.substr(skipped_read_length);
        this->CigarData.insert(this->CigarData.end(),
                                new_cigar_data.begin(), new_cigar_data.end());
        // merge consective ops
        for (i = 0, j = 0; j< this->CigarData.size();) {
            j ++;
            while(j != this->CigarData.size()) {
                if (this->CigarData[i].Type == this->CigarData[j].Type) {
                    this->CigarData[i].Length += this->CigarData[j].Length;
                    j++;
                }
                else {
                    this->CigarData[i+1] = this->CigarData[j ];
                    break;
                }
            }
            i ++;
        }
        this->CigarData.resize(i);
    }
    else {
        // Good news, first read overlapp second read
        // and simply do nothing
        ;
    }
    AlignmentExpertStarter aes;
    this->investigate(&aes);

    // DEBUG: output merged ones
    // fprintf(stdout, "Merged: \n");
    //output_bamalignment(first);
}

void JewelerAlignment::investigate(AlignmentExpert *ae) {
    // must perform the compatible test first!
    // justify whether the sequences contains the JewelerAlignment

    // not sure this is the case, but the coordinates are messed up
    // sometimes. TODO: read the samtools's specification, and make
    // the coordinates right.
    ae->initialize(this);

    int genome_start=this->get_start_position();
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


int JewelerAlignment::get_skipped_region(int skipped_genome_length,
                                         vector<CigarOp> &cigar_data,
                                         int &skipped_length) {

	vector<CigarOp>::const_iterator cigar_iter = this->CigarData.begin();
	vector<CigarOp>::const_iterator cigar_end  = this->CigarData.end();
	skipped_length = 0;
    cigar_data.clear();

	for (; cigar_iter != cigar_end; ++cigar_iter) {
		const CigarOp& op = (*cigar_iter);

		if (skipped_genome_length <= 0) {
			// the overlapped region has been processed.
			// directly add next op
			cigar_data.push_back(op);
			continue;
		}
		CigarOp new_op = (op);

		switch (op.Type) {

		case (Constants::BAM_CIGAR_MATCH_CHAR):
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR):
		case (Constants::BAM_CIGAR_MISMATCH_CHAR):
			// skip the overlapped region
			skipped_length += min(skipped_genome_length, (int) op.Length);
			// fall through
		case (Constants::BAM_CIGAR_REFSKIP_CHAR):
		case (Constants::BAM_CIGAR_DEL_CHAR):
			skipped_genome_length = skipped_genome_length - op.Length;
			new_op.Length = - skipped_genome_length;
			//fall trough

		case (Constants::BAM_CIGAR_PAD_CHAR):
			// not completely skipped, go to next opk
			if (skipped_genome_length >= 0) {
				continue;
			}
			cigar_data.push_back(new_op);
			break;

		case (Constants::BAM_CIGAR_INS_CHAR):
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR):
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			skipped_length += (int) op.Length;
			// simply do nothing, because these chars are from the
			// read, not the genome, so won't affect the skipped
			// region;
			break;
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			fprintf(stdout, "%s\n", message.c_str());
			exit(0);
		}
	}
	return 0;
}

