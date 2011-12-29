#include "alignment_glue.hpp"


int output_bamalignment(BamAlignment *al){
	fprintf(stdout,"%s\n",al->QueryBases.c_str());
	fprintf(stdout,"%d\n",al->Position);
	fprintf(stdout,"%s\n",get_cigar_string((*al)).c_str());
	return 0;
}
int AlignmentGlue::get_skipped_region(BamAlignment *al, int skipped_alignment_length, 
									  vector<CigarOp> &cigar_data, int &skipped_length){

	vector<CigarOp>::const_iterator cigar_iter = al->CigarData.begin();
	vector<CigarOp>::const_iterator cigar_end  = al->CigarData.end();
	skipped_length = 0;

	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);

		if (skipped_alignment_length <= 0) {
			// the overlapped region has been processed. 
			// directly add next op
			cigar_data.push_back(op);
			continue;
		}
		CigarOp new_op = (op);
		
		switch ( op.Type ) {

		case ( Constants::BAM_CIGAR_MATCH_CHAR ):
		case ( Constants::BAM_CIGAR_SEQMATCH_CHAR ): 
		case ( Constants::BAM_CIGAR_MISMATCH_CHAR ):
			// skip the overlapped region
			skipped_length += min(skipped_alignment_length, ( int ) op.Length);
			skipped_alignment_length = skipped_alignment_length - op.Length;
			new_op.Length = - skipped_alignment_length;

			//fall trough

		case ( Constants::BAM_CIGAR_DEL_CHAR):
		case ( Constants::BAM_CIGAR_PAD_CHAR):
		case ( Constants::BAM_CIGAR_REFSKIP_CHAR):

			// not completely skipped, go to next opk
			if (skipped_alignment_length >= 0 ){
				continue;
			}
			cigar_data.push_back(new_op);
			break;

		case ( Constants::BAM_CIGAR_INS_CHAR ):
		case ( Constants::BAM_CIGAR_SOFTCLIP_CHAR ):
		case ( Constants::BAM_CIGAR_HARDCLIP_CHAR) :
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


int AlignmentGlue::glue_paired_alignments(BamAlignment *first, BamAlignment *second){
	// DEBUG:
	// output_bamalignment(first);
	// output_bamalignment(second);
	int i;
	int first_start	 = first->Position;
	int first_end	 = first_start + first->AlignedBases.length();
	int second_start = second->Position;
	int second_end	 = second_start + second->AlignedBases.length();
	int overlapped	 = second_start - first_end;

	if (overlapped>=0){
		// Great news! no overlapped region between first and second
		// reads.
		first->Length += second->Length;

		first->Qualities += second->Qualities;
		string pad;
		for (i = 0; i < overlapped; i++){
			pad += 'N';
		}
		first->QueryBases += second->QueryBases;
		first->AlignedBases += pad + second->AlignedBases;
		if (overlapped > 0){
			first->CigarData.push_back(CigarOp('J',overlapped));
		}
		first->CigarData.insert(first->CigarData.end(), 
								second->CigarData.begin(), second->CigarData.end());
	}
	else if (second_end > first_end) {
		// bad news! Overlapped region need to be removed from the
		// second read. 
		int skipped_read_length = 0 , skipped_alignment_length = -overlapped;
		vector<CigarOp> new_cigar_data;
		get_skipped_region(second, skipped_alignment_length, 
						   new_cigar_data, skipped_read_length);
		first->Length += second->Length - skipped_read_length;
		first->QueryBases += second->QueryBases.substr(skipped_read_length);
		first->Qualities +=second->Qualities.substr(skipped_read_length);
		first->AlignedBases +=  second->AlignedBases.substr(skipped_read_length);
		first->CigarData.insert(first->CigarData.end(), 
								new_cigar_data.begin(), new_cigar_data.end());
		
	}
	else {
		// Good news, first read overlapp second read
		// and simply do nothing
		;
	}
	
	// DEBUG: output merged ones
	// fprintf(stdout, "Merged: \n");
	// output_bamalignment(first);
}

int AlignmentGlue::glue(vector<BamAlignment *> &in_reads, 
						vector<BamAlignment *> &new_reads,
						vector<BamAlignment *> &noused){
	int i;
	
	name2reads.clear();
	
	for ( i = 0; i < in_reads.size(); i++) {
		name2reads[ in_reads[i]->Name ].push_back(in_reads[ i ]);
	}

	for ( i = 0; i < in_reads.size(); i++) {
		if ( in_reads[i]->IsFirstMate() ) {
			if ( name2reads.find(in_reads[i]->Name) != name2reads.end()) {
				vector<BamAlignment *> alignments = name2reads[in_reads[i]->Name];
				if ( alignments.size() == 1 ) {
					// not paired, but still counted as a valid alignment.
					continue;
				}
				else if (alignments.size()==2){
					
					if (alignments[0]->IsFirstMate() && alignments[1]->IsSecondMate()) {
						glue_paired_alignments(alignments[0], alignments[1]);
						noused.push_back(alignments[1]);
						new_reads.push_back(alignments[0]);
					}
					else if (alignments[1]->IsFirstMate() && alignments[0]->IsSecondMate()) {
						glue_paired_alignments(alignments[1], alignments[0]);
						noused.push_back(alignments[0]);
						new_reads.push_back(alignments[1]);
					}
					else {
						//fprintf(stdout, "not properly mapped\n");
					}
				}
				else {
					//fprintf(stdout, 
					//"Oops, there are mulitple alignments within a gene?!\n");
					//for (int i = 0; i < alignments.size(); i++){
					//output_bamalignment(alignments[i]);
					//}
				}
			}
		}
	}
	
	return 0;
}
