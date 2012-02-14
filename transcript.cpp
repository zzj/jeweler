
#include "transcript.hpp"

Transcript::Transcript(){
	is_initialized=false;
}

int Transcript::load_seq(FastaReference * fr){
	seq = "";
	for ( int i=0; i < exon_start.size(); i++){
		seq+=fr->getSubSequence(chr, exon_start[i] - 1, exon_end[i]-exon_start[i]+1);
	}
	return 0;
}

int Transcript::insert_reads(JewelerAlignment *al ){
	reads.insert(al);
	return 0;
}

bool Transcript::is_aligned(JewelerAlignment *al ){
	return reads.find(al)!=reads.end();
}

int Transcript::get_next_exon(int start_pos, int start_seg = 0, int tolerate = 0) {
	while(!(exon_start[start_seg] - tolerate <= start_pos &&
			exon_end[start_seg] + tolerate >= start_pos)){
		start_seg++;
		if (start_seg>=exon_start.size())
			return NOT_FOUND;
	}
	return start_seg;
}

int Transcript::get_overlapped_alignment(JewelerAlignment *al , int &penalty, bool is_to_fix){

	penalty = 0;

	if (exon_start.size()!=exon_end.size()){
		fprintf(stderr,"the sizes does not match.\n");
		exit(0);
	}
	// Get start segments
	int start_seg = 0;
	int start_pos = al->Position + 1;

	int new_length = 0;

	int alignment_start = 0;
	string new_querybases;
	vector<CigarOp> new_cigar_data;
	int begin_seg, end_seg;
	int begin_err, end_err;
	int temp_length;
	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - ;
			// check wether the matched string belong to the same exon

		case ( Constants::BAM_CIGAR_MATCH_CHAR )    :
		case ( Constants::BAM_CIGAR_SEQMATCH_CHAR ) :
		case ( Constants::BAM_CIGAR_MISMATCH_CHAR ) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. 


			// two cases will be completely ignored
			// Case ONE: beginning and end of the matched region does
			// not belong to any exon region
			// Case TWO: beginning and end of the matched region
			// belong to two different exon region
			begin_seg = get_next_exon(start_pos);
			end_seg = get_next_exon(start_pos + op.Length -1);

			// testing
			// if ( al->Name == "UNC12-SN629_0154:8:2304:17730:106338#GGCTAC") {
			// 	fprintf(stdout, "%s\t%s\t%d\t%d\t%d\t%d\t%d\n", al->Name.c_str(),
			// 			get_cigar_string(*al).c_str(), al->Position + 1,
			// 			start_pos, begin_seg, end_seg,start_pos + op.Length -1);
			// 	fprintf(stdout, "%d\t%d\n",exon_start[begin_seg],exon_end[begin_seg]);
			// 	fprintf(stdout, "%d\t%d\n",exon_start[end_seg],exon_end[end_seg]);
			// }

			if  (( begin_seg == NOT_FOUND && begin_seg == end_seg)
				 ||( begin_seg != NOT_FOUND && end_seg != NOT_FOUND &&
					 begin_seg != end_seg) ) {

				// TODO: 
				// An anti-example the matched region can be covered
				// one or several exon regions.
				penalty += op.Length;
				start_pos += op.Length;
				new_cigar_data.push_back(CigarOp(Constants::BAM_CIGAR_REFSKIP_CHAR, 
												 op.Length));
			}
			else { 
				begin_err = 0; end_err = 0;

				if ( begin_seg == NOT_FOUND && end_seg != NOT_FOUND ){
					// the begin of the read exceeds the exon region
					begin_err = exon_start[end_seg] - start_pos;
					new_cigar_data.push_back(CigarOp(Constants::BAM_CIGAR_REFSKIP_CHAR, 
													 begin_err));
				}
				if ( begin_seg != NOT_FOUND && end_seg == NOT_FOUND ){
					// the end of the read exceeds the exon region
					end_err =  start_pos + op.Length -1 - exon_end[ begin_seg ];
				}
				start_pos += op.Length  ;

				temp_length = op.Length - begin_err - end_err;

				new_cigar_data.push_back(CigarOp(op.Type, temp_length));

				new_querybases += al->QueryBases.substr(alignment_start + begin_err,
														temp_length);
				new_length += temp_length;
				alignment_start += op.Length;
				if (end_err > 0 ){
					new_cigar_data.push_back(CigarOp(Constants::BAM_CIGAR_REFSKIP_CHAR, 
													 end_err));
				}
				penalty += begin_err + end_err;
			}
					 
			break;
			
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			if (get_next_exon(start_pos) != NOT_FOUND) {
				new_cigar_data.push_back(op);
				new_length += op.Length;
				new_querybases += al->QueryBases.substr(alignment_start, op.Length);
				alignment_start += op.Length;
			}
		    else {
				penalty += op.Length;
			}
			break;
			
		// for 'N', goto next exon
		// Only 'N' should be presented as an intron, because 'D' and 'P' are not
		// defined in RNA-seq data based on the Samtools document.
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case ('J'):
			start_pos += op.Length;
			new_cigar_data.push_back(op);
			break;
			// invalid CIGAR op-code

			
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			fprintf(stderr, "%s\n", message.c_str());
			exit(0);

		}
	}
	
	if (penalty > 0 && is_to_fix) {
		// fprintf(stdout, "Find a read need to be fixed!\n");
		// output_segments();
		// fprintf(stdout,"%s\n",al->Name.c_str());
		// fprintf(stdout,"%s\n",al->QueryBases.c_str());
		// fprintf(stdout,"%d\n",al->Position + 1);
		// fprintf(stdout,"%s\n",get_cigar_string((*al)).c_str());
		al->Length = new_length;
		al->CigarData = new_cigar_data;
		al->QueryBases = new_querybases;

		cigar_trim(al);

		// fprintf(stdout,"%s\n",get_cigar_string((*al)).c_str());
		// fprintf(stdout,"%s\n",al->QueryBases.c_str());
		// fprintf(stdout,"%d\n",al->Position + 1);


	}

	return true;
}

bool Transcript::is_compatible(JewelerAlignment *al , int tolerate , bool debug){
	// justify whether the sequences contains the JewelerAlignment
	// if it contains most of the read, this will adjust the read to
	// make it completely compatible 
	// check it out with Transcript::tolerate

	if (exon_start.size()!=exon_end.size()){
		fprintf(stderr,"the sizes does not match.\n");
		exit(0);
	}
	// Get start segments
	int start_seg=0;
	int start_pos=al->Position + 1;

	int err;
	if ( (start_seg = get_next_exon (  start_pos, start_seg) ) == NOT_FOUND ){
		return false;
	}
	
	int new_length = 0;
	string new_querybases = "";

	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - ;
			// check wether the matched string belong to the same exon

		case ( Constants::BAM_CIGAR_MATCH_CHAR )    :
		case ( Constants::BAM_CIGAR_SEQMATCH_CHAR ) :
		case ( Constants::BAM_CIGAR_MISMATCH_CHAR ) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. 
			start_pos += op.Length;
			if ( !( exon_start[start_seg] - tolerate <= start_pos &&
					exon_end[start_seg] + tolerate >= start_pos -1 )
				 ) {// the end of last alignment 
				if (debug){
					fprintf(stdout, "not at the same exon for matching region, current position = %d\n",
							start_pos);
				}

			    return false;
			}
			break;
			
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			break;

		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
		// for 'N', goto next exon
		// Only 'N' should be presented as an intron, because 'D' and 'P' are not
		// defined in RNA-seq data based on the Samtools document.
			err = (start_pos - 1 - exon_end[ start_seg ]);
			if (  abs(err) <= tolerate){
				// start postion at the end of current exon
				start_pos += op.Length;
				start_seg ++;
				err = start_pos - exon_start[ start_seg ];
				if ( abs(err) > tolerate) {
					// not at the beginning of the exon
				if (debug){
					fprintf(stdout, "not at start of next exon, current position = %d\n",
							start_pos);
				}

					return false;
				}
			}
			else {
				// not at the end of exon
				if (debug){
					fprintf(stdout, "not at end of current exon, current position = %d\n",
							start_pos);
				}
				return false;
			}
			break;
			
		case (Constants::BAM_CIGAR_DEL_CHAR) :
			start_pos += op.Length;
			break;
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case ('J') :
			start_pos += op.Length;
			if ( (start_seg = get_next_exon (  start_pos, start_seg) ) == NOT_FOUND ){
				if (debug){
					fprintf(stdout, "did not find next exon, current position = %d\n",
							start_pos);
				}
				return false;
			}
			break;
			// invalid CIGAR op-code
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			return false;
		}
	}
	
	if (debug){
		fprintf(stdout, "everything is fine\n");
	}
	return true;
}

int Transcript::output_segments(){
	for (int i=0;i<exon_start.size();i++){
		fprintf(stdout,"(%d,%d)",exon_start[i],exon_end[i]);
	}
	fprintf(stdout,"\n");
}


string Transcript::get_query_aligned_seq(JewelerAlignment * al){
	// must perform the compatible test first!
	// justify whether the sequences contains the JewelerAlignment

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
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
			// insertion, soft clipnot aligned 
			start_pos+=op.Length;
			break;


		// for 'D' - write gap character
		// for 'N' - write N's, skip bases in original query sequence
		// for 'P' - write padding character			
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
		case ('J'):
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

int Transcript::match_alleles(JewelerAlignment *al, int &total_alleles, 
							  vector<int> &transcript_aligned_locations,
							  vector<int> &read_aligned_locations,
							  vector<int> &matched_alleles, 
							  vector<int> &mismatches,
							  vector<char> & read_mismatch_qualities,
							  vector<char> & mismatchars){
	int i;
	string transcript_seq = get_transcript_aligned_info<string>(al, get_seq_info);
	transcript_aligned_locations =
		get_transcript_aligned_info<vector<int> >(al, get_transcript_location_info);
	read_aligned_locations =
		get_transcript_aligned_info<vector<int> >(al, get_read_location_info);
	string query_seq=get_query_aligned_seq(al);

	alleles.clear();

	matched_alleles.clear();
	mismatches.clear();
	mismatchars.clear();
	total_alleles=0;

	if (transcript_seq.size()!=query_seq.size()){
		
		fprintf(stderr,"Transcript: %d\t%s\n", transcript_seq.size(), transcript_seq.c_str());
		fprintf(stderr,"Read:     : %d\t%s\n", query_seq.size(), query_seq.c_str());
		fprintf(stderr,"%s\t%d\t%s\n", 
				al->Name.c_str(), 
				al->Position + 1, 
				get_cigar_string(*al).c_str());

		output_segments();
		
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
			read_mismatch_qualities.push_back(al->Qualities[i]);
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
		if (genome_location >= exon_start[i]-tolerate  && 
			genome_location <= exon_end[i] + tolerate){
			return i;
		}
	}

	return NOT_FOUND;
}

int Transcript::register_allele_read(JewelerAlignment *al){
	int total_alleles;
	vector<int> matched_alleles;
	vector<int> locations;
	vector<int> read_locations;
	vector<int> mismatches;
	vector<char> read_mismatches;
	vector<int> matched_exons;
	vector<char> mismatchars;
	int num_matched_alleles;
	int i;
	match_alleles(al,total_alleles,locations, read_locations, matched_alleles,
				  mismatches,read_mismatches,mismatchars);
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

int Transcript::register_read(JewelerAlignment *al){

	vector<int> matched_exons;
	int i;
	matched_exons=get_transcript_aligned_info<vector<int> > (al,get_exon_info);
	for (i = 0; i < matched_exons.size(); i++){
		reads_per_exon[matched_exons[i]].insert(al);
	}
	reads.insert(al);
	
	return 0;
}


int get_seq_info(Transcript * ti, JewelerAlignment *al,
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
int get_transcript_location_info(Transcript * ti, JewelerAlignment *al,
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
int get_read_location_info(Transcript * ti, JewelerAlignment *al,
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
		ret.push_back(alignment_start+i);
	}
	return 0;
}

int get_exon_info(Transcript * ti,  JewelerAlignment * al, 
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

int insert_mismatch_info(Transcript *ti,  JewelerAlignment * al, 
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
