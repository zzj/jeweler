
#include "transcript.hpp"

Transcript::Transcript(){
	is_initialized=false;
}


int Transcript::insert_aligned_reads(BamAlignment *al ){
	aligned_reads.insert(al);
	return 0;
}

bool Transcript::is_aligned(BamAlignment *al ){
	return aligned_reads.find(al)!=aligned_reads.end();
}


bool Transcript::is_compatible(BamAlignment *al ){
	// justify whether the sequences contains the BamAlignment

	if (exon_start.size()!=exon_end.size()){
		fprintf(stderr,"the sizes does not match.\n");
		exit(0);
	}
	// Get start segments
	int start_seg=0;
	int start_pos=al->Position+1;
	while(!(exon_start[start_seg]<=start_pos &&
			exon_end[start_seg]>=start_pos)){
		start_seg++;
		if (start_seg>=exon_start.size())
			return false;
	}

	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - write bases
		case (Constants::BAM_CIGAR_MATCH_CHAR)    :
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
		case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. 
			start_pos+=op.Length;
			if (!(exon_start[start_seg]<=start_pos &&
				  exon_end[start_seg]>=start_pos-1)) // the end of
													 // last alignment 
			    return false;
			break;
			
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
			// fall through
			break;
		// for 'S' - soft clip, do not write bases
		// but increment placeholder 'k'
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
			break;
			
		// for 'D' - write gap character
		// for 'N' - write N's, skip bases in original query sequence
		// for 'P' - write padding character			
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
			

			// only increase one segment if it matches the end of exon
		    // otherwise, does not increase segments, because it may
		    // be deletion
			if (start_pos==exon_end[start_seg]+1){
				start_pos+=op.Length;
				start_seg++;
				if (!(exon_start[start_seg]==start_pos )){
					//fprintf(stdout,"no matched! %d\t %d\t\n", start_pos, exon_start[start_seg]);
					return false;
				}
				else {
					//fprintf(stdout,"matched! %d\t %d\t\n", start_pos, exon_start[start_seg]);
				}
			}
			else 			start_pos+=op.Length;
			break;
			
			// for 'H' - hard clip, do nothing to AlignedBases, move to next op
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			break;
			
			// invalid CIGAR op-code
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			return false;
		}
	}
	
	return true;
}






int Transcript::output_segments(){
	for (int i=0;i<exon_start.size();i++){
		fprintf(stdout,"(%d,%d)",exon_start[i],exon_end[i]);
	}
	fprintf(stdout,"\n");
}


string Transcript::get_transcript_aligned_seq(BamAlignment * al){
	// must perform the compatible test first!
	// justify whether the sequences contains the BamAlignment

	// not sure this is the case, but the coordinates are messed up
	// sometimes. TODO: read the samtools's specification, and make
	// the coordinates right.
	int start_pos=al->Position+1;
	int start;
	string aligned_string;
	
	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - write bases
		case (Constants::BAM_CIGAR_MATCH_CHAR)    :
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
		case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. 
			start=get_transcript_location(start_pos);
			if (start==-1){
				fprintf(stderr,"cannot find transcript postition, did you run Transcript::is_compatible ?\n %d\n%d\n", start_pos,genome_pos.size());

				output_segments();
				exit(0);
			}
			aligned_string+=seq.substr(start,op.Length);
			start_pos+=op.Length;
			break;
			
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
				// fall through
			break;
		// for 'S' - soft clip, do not write bases
		// but increment placeholder 'k'
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :

			break;
			
		// for 'D' - write gap character
		// for 'N' - write N's, skip bases in original query sequence
		// for 'P' - write padding character			
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
			start_pos+=op.Length;
			break;
			
			// for 'H' - hard clip, do nothing to AlignedBases, move to next op
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			break;
			
			// invalid CIGAR op-code
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			return false;
		}
	}
	return aligned_string;
	
}

string Transcript::get_query_aligned_seq(BamAlignment * al){
	// must perform the compatible test first!
	// justify whether the sequences contains the BamAlignment

	// not sure this is the case, but the coordinates are messed up
	// sometimes. TODO: read the samtools's specification, and make
	// the coordinates right.
	int start_pos=0;
	int start;
	string aligned_string;
	string query_seq=al->QueryBases;
	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - write bases
		case (Constants::BAM_CIGAR_MATCH_CHAR)    :
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
		case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. 
			aligned_string+=query_seq.substr(start_pos,op.Length);
			start_pos+=op.Length;
			break;
			
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
			// insertion, not aligned 
			start_pos+=op.Length;
			break;
		// for 'S' - soft clip, do not write bases
		// but increment placeholder 'k'
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
			break;
			
		// for 'D' - write gap character
		// for 'N' - write N's, skip bases in original query sequence
		// for 'P' - write padding character			
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
			break;
			
			// for 'H' - hard clip, do nothing to AlignedBases, move to next op
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			break;
			
			// invalid CIGAR op-code
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			return false;
		}
	}
	return aligned_string;
	
}

int Transcript::match_alleles(BamAlignment *al, int &total_alleles, int &num_matched_alleles,
							  vector<int> &matched_alleles){
	int i;
	string transcript_seq=get_transcript_aligned_seq(al);
	string query_seq=get_query_aligned_seq(al);
	int transcript_location=get_transcript_location(al->Position+1);
	alleles.clear();
	total_alleles=0;
	num_matched_alleles=0;
	if (transcript_seq.size()!=query_seq.size()){
		fprintf(stderr,"The aligned sequences from the transcript and the query do not match\n");
		exit(0);
	}
	// TODO: increase the performance
	for (i=0;i<transcript_seq.size();i++){
		// is it a SNP?
		if (is_allele(transcript_location+i)){
			total_alleles++;

			if (transcript_seq[i]!=query_seq[i]){
				;
			}
			else {
				num_matched_alleles++;
				matched_alleles.push_back(transcript_location+i);
			}
		}
	}
	return 0;
}



bool Transcript::is_allele(int transcript_location){
	return find(snp_pos.begin(),snp_pos.end(),transcript_location) !=snp_pos.end();
}


char Transcript::get_allele_char(int transcript_location){
	return alleles[lower_bound(snp_pos.begin(),snp_pos.end(),transcript_location) 
				   -snp_pos.begin()];
}
int Transcript::get_allele_exon(int transcript_location){
	return allele_exon[lower_bound(snp_pos.begin(),snp_pos.end(),transcript_location) 
				   -snp_pos.begin()];
}


int Transcript::get_transcript_location(int genome_location){

	int ret = lower_bound(genome_pos.begin(),genome_pos.end(),genome_location) - genome_pos.begin();
	if (ret==genome_pos.size()){
		return -1;
	}
	return ret;
}


int Transcript::register_read(BamAlignment *al){
	int total_alleles;
	vector<int> matched_alleles;
	int num_matched_alleles;
	int i;
	match_alleles(al,total_alleles,num_matched_alleles,matched_alleles);
	for (i=0;i<num_matched_alleles;i++){
		num_info_reads_per_exon[get_allele_exon(matched_alleles[i])]++;
	}

	aligned_reads.insert(al);

	if (num_matched_alleles>0)
		allele_aligned_reads.insert(al);
	
	return 0;
}
