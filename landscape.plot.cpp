

#include "landscape.plot.hpp"

LandscapePlot::LandscapePlot(Transcript * maternal_transcript,
							 Transcript * paternal_transcript,
							 set<BamAlignment *> & unknown_reads){
	
	int i,j;

	if (maternal_transcript->seq.size()!=paternal_transcript->seq.size()){
		fprintf(stderr,"Maternal and paternal transcripts must be the same one");
		exit(0);
	}
	transcript_name=maternal_transcript->name;
	int size=maternal_transcript->seq.size();
	
	unknown.resize(size,0);
	paternal.resize(size,0);
	maternal.resize(size,0);
	is_snp.resize(size,0);
	exon_jump.resize(size,0); // the last position that an exon ends
	this->num_maternal=maternal_transcript->allele_aligned_reads.size();
	this->num_paternal=paternal_transcript->allele_aligned_reads.size();
	this->num_unknown=unknown_reads.size();

	for (i=0;i<maternal_transcript->snp_pos.size();i++){
		is_snp[maternal_transcript->snp_pos[i]]=1;
	}
	int current_pos=0;
	for (j=0;j<maternal_transcript->exon_start.size();j++){
		//exon_end is inclusive
		current_pos=current_pos+maternal_transcript->exon_end[j]-maternal_transcript->exon_start[j]+1;
		exon_jump[current_pos-1]=1;
	}

	add_transcript_to_landscape(maternal_transcript, 
								maternal_transcript->allele_aligned_reads,
								maternal);
	add_transcript_to_landscape(paternal_transcript, 
								paternal_transcript->allele_aligned_reads,
								paternal);

	add_transcript_to_landscape(paternal_transcript /* either one is
													   fine */,  
								unknown_reads,
								unknown);


}

int LandscapePlot::add_transcript_to_landscape(Transcript * transcript,
										  set<BamAlignment *> &reads,
										  vector<int> & coverage){
	for (auto iter=reads.begin();
		 iter!=reads.end();
		 iter++){
		add_coverage(transcript, *iter, coverage);
	}
	return 0;
}


int LandscapePlot::add_coverage(Transcript * transcript, BamAlignment *al, vector<int> &coverage){
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
	int i;
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - write bases
		case (Constants::BAM_CIGAR_MATCH_CHAR)    :
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
		case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. 
			start=transcript->get_transcript_location(start_pos);
			for (i=start;i<start+op.Length;i++){
				coverage[i]++;
			}
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
	return 0;
}

int LandscapePlot::generate_landscape_plot(FILE * info, FILE * output){
	int i;
	fprintf(info,"%s\t%d\t%d\t%d\n",transcript_name.c_str(),num_unknown,num_paternal,num_maternal);
	
	for (i=0;i<unknown.size();i++){
		fprintf(output,"%d\t%d\t%d\t%d\t%d\n", unknown[i], paternal[i], maternal[i], is_snp[i], exon_jump[i]);
	}
			
}
