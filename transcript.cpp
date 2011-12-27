
#include "transcript.hpp"

Transcript::Transcript(){
	is_initialized=false;
}

int Transcript::insert_reads(BamAlignment *al ){
	reads.insert(al);
	return 0;
}

bool Transcript::is_aligned(BamAlignment *al ){
	return reads.find(al)!=reads.end();
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

int Transcript::match_alleles(BamAlignment *al, int &total_alleles, 
							  vector<int> &transcript_aligned_locations,
							  vector<int> &matched_alleles, 
							  vector<int> &mismatches,
							  vector<char> & mismatchars){
	int i;
	string transcript_seq = get_transcript_aligned_info<string>(al, get_seq_info);
	transcript_aligned_locations =
		get_transcript_aligned_info<vector<int> >(al, get_location_info);
	string query_seq=get_query_aligned_seq(al);
	
	alleles.clear();

	matched_alleles.clear();
	mismatches.clear();
	mismatchars.clear();
	total_alleles=0;

	if (transcript_seq.size()!=query_seq.size()){
		fprintf(stderr,"The aligned sequences from the transcript and the query do not match\n");
		exit(0);
	}
	
	// TODO: increase the performance
	for (i=0;i<transcript_seq.size();i++){
		// is it a SNP?
		if (is_allele(transcript_aligned_locations[i])){
			total_alleles++;
			if (transcript_seq[i]==query_seq[i]){
				matched_alleles.push_back(transcript_aligned_locations[i]);
				continue;
			}
		}
		// is a mismatch			
		if (transcript_seq[i] != query_seq[i]){
			mismatches.push_back(transcript_aligned_locations[i]);
			mismatchars.push_back(query_seq[i]);
		}
	}
	if (total_alleles>0 && matched_alleles.size()==0){
		//fprintf(stdout, "%d\t%d\n", total_alleles, num_matched_alleles);
		//fprintf(stdout, "transcript_seq: %s\n", transcript_seq.c_str());
		//fprintf(stdout, "query_seq:      %s\n", query_seq.c_str());
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

int Transcript::get_transcript_exon(int genome_location){
	
	int i;
	for (i=0;i<exon_start.size();i++){
		if (genome_location>=exon_start[i] && 
			genome_location<=exon_end[i]){
			return i;
		}
	}

	return NOT_FOUND;
}

int Transcript::register_allele_read(BamAlignment *al){
	int total_alleles;
	vector<int> matched_alleles;
	vector<int> locations;
	vector<int> mismatches;
	vector<int> matched_exons;
	vector<char> mismatchars;
	int num_matched_alleles;
	int i;
	match_alleles(al,total_alleles,locations, matched_alleles,mismatches,mismatchars);
	num_matched_alleles=matched_alleles.size();
	if (num_matched_alleles>0){
		allele_reads.insert(al);
		for (i=0;i<num_matched_alleles;i++){
			int exon_id=get_allele_exon(matched_alleles[i]);
			num_info_reads_per_exon[exon_id]++;
			allele_reads_per_exon[exon_id].insert(al);
		}
		
	}
	register_read(al);
	return 0;
}

int Transcript::register_read(BamAlignment *al){

	vector<int> matched_exons;
	int i;
	matched_exons=get_transcript_aligned_info<vector<int> > (al,get_exon_info);
	for (i = 0; i < matched_exons.size(); i++){
		reads_per_exon[matched_exons[i]].insert(al);
	}
	reads.insert(al);
	
	
	return 0;
}


int get_seq_info(Transcript * ti, BamAlignment *al,
				 int genome_start, int alignment_start, int length, 
				 string &ret){
	int transcript_start=ti->get_transcript_location(genome_start);
	if (transcript_start==NOT_FOUND){
		fprintf(stderr,"cannot find transcript postition, did you run Transcript::is_compatible ?\n %d\n%d\n", genome_start,ti->genome_pos.size());
		ti->output_segments();
		exit(0);
	}
	ret+=ti->seq.substr(transcript_start,length);
	return 0;
}
int get_location_info(Transcript * ti, BamAlignment *al,
					  int genome_start, int alignment_start, int length, 
					  vector<int>  &ret){
	int transcript_start = ti->get_transcript_location(genome_start);
	int i;
	if (transcript_start == NOT_FOUND){
		fprintf(stderr,"cannot find transcript postition, did you run Transcript::is_compatible ?\n %d\n%d\n", genome_start,ti->genome_pos.size());
		ti->output_segments();
		exit(0);
	}
	for (i = 0 ; i < length ; i++){
		ret.push_back(transcript_start+i);
	}
	return 0;
}

int get_exon_info(Transcript * ti,  BamAlignment * al, 
				  int genome_start, int alignment_start, int length, 
				  vector<int> &ret){

	int transcript_start=ti->get_transcript_location(genome_start);
	if (transcript_start==NOT_FOUND){
		fprintf(stderr,"cannot find transcript postition, did you run Transcript::is_compatible ?\n %d\n%d\n", genome_start,ti->genome_pos.size());
		
		ti->output_segments();
		exit(0);
	}
	ret.push_back(ti->get_transcript_exon(genome_start));
	return 0;
}

int insert_mismatch_info(Transcript *ti,  BamAlignment * al, 
						 int genome_start, int alignment_start, int length, 
						 vector<int> &mismatches){
	int i;
	int transcript_start=ti->get_transcript_location(genome_start);
	if (transcript_start==NOT_FOUND){
		fprintf(stderr,"cannot find transcript postition, did you run Transcript::is_compatible ?\n %d\n%d\n", genome_start,ti->genome_pos.size());
		
		ti->output_segments();
		exit(0);
	}
	string genome_seq    = ti->seq.substr(transcript_start, length);
	string alignment_seq = al->QueryBases.substr(alignment_start, length);
	for ( i = 0; i < genome_seq.size(); i++){
		if (genome_seq[i] != alignment_seq[i]){
			mismatches.push_back( genome_start + i );
		}
	}
	return 0;
}

int Transcript::add_transcript_to_graph(Graph *graph, vector<Path> &records){
	vector<ExonNode *> exon_chain;
	bool is_valid=true;
	exon_chain.resize(exon_start.size(),NULL);
	
	for (int i=0;i<exon_start.size();i++){
		int exon_origin=EXON_NO_INFO;
		if (num_alleles_per_exon[i]>0){
			if (num_info_reads_per_exon[i]>0){
				exon_origin=origin;
			}
			else {
				// do not insert this exon, because there is no reads
				// found to be informative, though there are several
				// alleles
				is_valid=false;
				continue;
			}
		}
		exon_chain[i]=graph->add_exon_node(exon_start[i],exon_end[i],
										   exon_origin, reads_per_exon[i],
										   allele_reads_per_exon[i]);
		if (i==0) {
			continue;
		}
		// TODO: fix the number of reads
		if (exon_chain[i-1]!=NULL){
			graph->add_edge(exon_chain[i-1],exon_chain[i],1);
		}
	}

	records.push_back(Path(exon_chain));
	return 0;
}

bool Transcript::is_equal(Transcript *t){
	if (t->exon_start.size()!=this->exon_start.size()){
		return false;
	}
	for (int i=0;i<t->exon_start.size();i++){
		if (t->exon_start[i]!=this->exon_start[i]){
			return false;
		}
		if (t->exon_end[i]!=this->exon_end[i]){
			return false;
		}
	}
	return true;
}
