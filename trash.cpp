	// record all mismatches in transcript
	// this function will change transcript. Please be careful about it.
	// make sure for each read, it is be called once. Otherwise,
	// the statistics for mismaches are wrong. 
	int annotate_mismatch_pos(Transcript * t, rna_read_query * rrq);


int jeweler::annotate_mismatch_pos(Transcript *t, rna_read_query * rrq)
{
	int num_mismatches = 0;
	int i,j;
	for( i = 0; i < rrq->query_start.size(); ++i)
	{
		int qstr = rrq->query_start[i];
		int tstr = rrq->target_start[i];
		int len  = rrq->block_size[i];
		
		for( j = 0; j < len; ++j)
		{

			if(t->seq[tstr+j]!=rrq->seq[qstr+j]) // find a mismatch or an allele
			{
									
				num_mismatches++;

				if(find(t->snp_pos.begin(), t->snp_pos.end(), tstr+j)!=t->snp_pos.end())
					continue;
				else
				{
					//printf("%s %s %d %d %d\n",t->name.c_str(),rrq->name.c_str(),tstr+len,t->noninformative_mismatches.size(),t->seq.size());
					if (tstr+j>t->noninformative_mismatches.size()) {
						fprintf(stderr,"out of range!\n");
					}
					t->noninformative_mismatches[tstr+j] += 1;
					// revised by Weibo
					switch(t->seq[tstr+j])
					{
					case 'A':
					case 'a':
						t->Anoninformative_mismatches[tstr+j] += 1;
						break;
					case 'C':
					case 'c':
						t->Cnoninformative_mismatches[tstr+j] += 1;
						break;
					case 'G':
					case 'g':
						t->Gnoninformative_mismatches[tstr+j] += 1;
						break;
					case 'T':
					case 't':
						t->Tnoninformative_mismatches[tstr+j] += 1;
					}
				}
			}
		}
	}
	//if (num_mismatches < rrq->mismatch)
	if (false)
	{

	    // the reason this happened is due to the randomness of blat
		// sometimes, blat does not return the optimal solution
		// so when it aligned to the first transcript, it may not 
		// align both pair properly.
		// however, when it align to the second transcript, it may align
		// both properly.
		// then it may cause this warning. 
		// the merged read can be aligned better to the improper aligned transcript, 
		// which has been discarded before.  
		fprintf(stderr, "%s\t%s\n", t->name.c_str(),t->seq.c_str());
		fprintf(stderr, "%s\t%s\n", rrq->name.c_str(),rrq->seq.c_str());
		fprintf(stderr, "%d %d %d %d\n",
				num_mismatches, rrq->mismatch,
				rrq->query_start.size(),count_mismatches(t,rrq));
		for (i=0;i<rrq->query_start.size();i++){
			fprintf(stderr,"%s\n",rrq->seq.substr(rrq->query_start[i],rrq->block_size[i]).c_str());
			fprintf(stderr,"%d\n",rrq->block_size[i]);
		}
		fprintf(stderr, "WARNING: inconsistent mismatches at query %s \n", rrq->name.c_str());
	}

	return num_mismatches;

}


int jeweler::generate_landscape(TranscriptInfo * ti,
								vector<Transcript *> &ref,
								map<rna_read_key,rna_read_query *> &queries){
	
	int i,j,k,m;
	string filename=string(ti->folder+ti->gene_id+"landscape.plot.info");
	FILE *fd=file_open(filename.c_str(),"w+");
	for ( i=0;i<ref.size();i++){
		int size=ref[i]->seq.size();
		int num_unknown=0,num_paternal=0,num_maternal=0;

		vector<int> unknown(size,0);
		vector<int> paternal(size,0);
		vector<int> maternal(size,0);
		vector<int> is_snp(size,0);
		// the last position that a  exon ends
		vector<int> exon_jump(size,0);
		
		for (j=0;j<ref[i]->snp_pos.size();j++) {
			is_snp[ref[i]->snp_pos[j]]=1;
		}

		int current_pos=0;
		for (j=0;j<ref[i]->exon_start.size();j++){
			//exon_end is inclusive
			current_pos=current_pos+ref[i]->exon_end[j]-ref[i]->exon_start[j]+1;
			exon_jump[current_pos-1]=1;
		}

		for (auto iter=queries.begin();iter!=queries.end();iter++){
			if (iter->first.transcript_name==ref[i]->name){

				vector<int> *target=NULL;

				if (iter->second->source_id==0){
					num_unknown++;
					target=&unknown;
				}
				if (iter->second->source_id==1){
					num_paternal++;
					target=&paternal;
				}
				if (iter->second->source_id==2){
					num_maternal++;
					target=&maternal;
				}
				if(target==NULL){
					fprintf(stderr,"No source id?!\n");
					exit(0);
				}

				for (k=0;k<iter->second->target_start.size();k++){
					for (m=0;m<iter->second->block_size[k];m++){
						if (iter->second->target_start[k]+m>target->size()){
							fprintf(stderr,"out of range\n");
							exit(0);
						}
						target->at(iter->second->target_start[k]+m)++;
					}
				}

				if (iter->second->first_end >0 && iter->second->first_end < iter->second->second_start)
				{
					for(int i = iter->second->first_end+1; i < iter->second->second_start; ++i)
					{
						target->at(i)++;
					}
				}
			}
		}

		for (k=0;k<ref[i]->snp_pos.size();k++){
			int target=ref[i]->snp_pos[k];
			//printf("%d\t%d\t%d\n",unknown[target],paternal[target],maternal[target]);
 		}
		
		if (i==0){
			string filename=string(ti->folder+ti->gene_id+"landscape.plot.info.meta");
			FILE *fd=file_open(filename.c_str(),"w+");
			fprintf(fd,"%d\t%d\t%d\n",num_unknown,num_paternal,num_maternal);
			fclose(fd);
		}


		for (k=0;k<ref[i]->seq.size();k++){
			// revised by weibo
			int target=k;
			fprintf(fd,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", 
				ref[i]->name.c_str(),target, unknown[target],paternal[target],maternal[target],
					is_snp[target],exon_jump[target],ref[i]->noninformative_mismatches[k],
					ref[i]->genome_pos[k], ref[i]->Anoninformative_mismatches[k],
					ref[i]->Cnoninformative_mismatches[k],
					ref[i]->Gnoninformative_mismatches[k],
					ref[i]->Tnoninformative_mismatches[k],
					ref[i]->chr.c_str());
		}
	}
	fclose(fd);
	return 0;
}
	int generate_landscape(TranscriptInfo * ti,
						   vector<Transcript *> &ref,
						   map<rna_read_key,rna_read_query *> &queries);


	// load read data,
	int load_read_data(TranscriptInfo * ti, 
					   vector<Transcript *> &ptrans,
					   vector<Transcript *> &mtrans,
					   map<rna_read_key,rna_read_query *>& queires);
	int merge_paired_reads(vector<Transcript *> &ref,
						  vector<rna_read_query *> & blat_result);
	//  [outdated] add queries, match the read sequences to the transcript's sequences. 
	// If the reads are not flipped and reversed, do that. 
	int add_queries(vector<Transcript *> &ref, 
					multimap<string,string> &srmap,
					vector<rna_read_query *> &pqueries,
					map<rna_read_key,rna_read_query *>& queires);

	// return wether the read map to the transcript. 
	bool match_snp(Transcript * t, rna_read_query * rrq);

	// 
	int label_mismatches_perbase(vector<Transcript *> &ptrans, 
								 vector<Transcript *> &mtrans,
								 map<rna_read_key,rna_read_query *>& queires);
	// given a set of transcripts, check whether the read can be aligned to
	// the transcript by matching all the SNPs in the transcript. 
	int identify_sources(vector<Transcript *> &source,
						 map<rna_read_key,rna_read_query *> &queries,
						 int source_id);
	// return the number of mismatches
	int count_mismatches(Transcript *t, rna_read_query *rrq);
	

int jeweler::count_mismatches(Transcript *t, rna_read_query  *rrq){
	if (t->name != rrq->target){
		fprintf(stderr,"transcript names do not match with each other\n");
		exit(0);
	}
	int i,j,mismatches=0;

	for (i=0;i<rrq->block_size.size();i++){
		for (j=0;j<rrq->block_size[i];j++){
			if (rrq->seq[rrq->query_start[i]+j] != t->seq[rrq->target_start[i]+j]){
				mismatches++;
			}
		}
	}
	return mismatches;
}

int jeweler::merge_paired_reads(vector<Transcript *> &ref,
								vector<rna_read_query *> &bam_result){
	
	int i,j,k;

	// first, map all reads with same read id and same target transcript id 
	// to the same vector
	
	map<rna_read_key, vector<rna_read_query *> > read_id2queries;
	map<rna_read_key, vector<rna_read_query *> > :: iterator iter;

	for (i=0;i<bam_result.size();i++){
		vector<string> strs;
		string read_id,flag_field;
		// first field is read_id
		// second field's first field is flag
		// second field's second field is rna_read_query
		if (bam_result[i]->is_ignored) {
			continue;
		}
		rna_read_key  rrk = rna_read_key(bam_result[i]->target,bam_result[i]->name);
		iter=read_id2queries.find(rrk);
		
		if (iter==read_id2queries.end()){
			vector<rna_read_query *> vrrqs;
			vrrqs.push_back(bam_result[i]);
			read_id2queries.insert(make_pair(rrk, vrrqs));
		}
		else {
			vector<rna_read_query *> &vrrqs=iter->second;
			bool has_conflict=false;
			for (j=0;j<vrrqs.size();j++){

				if (vrrqs[j]->flag_field!=bam_result[i]->flag_field){
					continue;
				}
				else {
					has_conflict=true;
					if (is_better_alignment(bam_result[i],vrrqs[j])){
						vrrqs[j]=bam_result[i];
					}
					// no matter this is a good alignment or not
					// break the loop
					break;
				}
			}
			if (!has_conflict){
				vrrqs.push_back(bam_result[i]);
			}
		}

	}

	// for all reads that can be mapped back to the same transcript and read
	// merge them, and generate new bam result
	bam_result.clear();
	
	for (iter=read_id2queries.begin();iter!=read_id2queries.end();iter++){
		vector<rna_read_query *>& vrrqs=iter->second;
		if (vrrqs.size()==1){
		}
		else if (vrrqs.size()==2) {

			if (vrrqs[0]->target_start>vrrqs[1]->target_start){
				swap(vrrqs[0],vrrqs[1]);
			}
			// using s simple rule. 
			// append the non-overlapping fields 
			// from the second read to the first one
			// TODO: hanle more complicated cases
			// the start position in the second read
			// that does not overlap with the first one
			
			
			// get the last mapped position in the target_transcript
			// from the first read
			// the block is left-closed, and right-open
			int l = vrrqs[0]->target_start.size()-1;
			int last_in_first_tran=vrrqs[0]->target_start[l]+vrrqs[0]->block_size[l];
			int last_in_first_read=vrrqs[0]->query_start[l]+vrrqs[0]->block_size[l];

			vrrqs[0]->first_end = last_in_first_tran;
			vrrqs[0]->second_start = vrrqs[1]->target_start[0];
			
			for (i=0;i<vrrqs[1]->target_start.size();i++){
				int tran_block_end=vrrqs[1]->target_start[i]+vrrqs[1]->block_size[i];
				if (tran_block_end<=last_in_first_tran){
					continue;
				}
				int tran_block_start=vrrqs[1]->target_start[i];
				// no overlap, happy, just append to the read
				int shift=0;
				// overlapped 
				if (tran_block_start<last_in_first_tran){
						shift=last_in_first_tran-tran_block_start;
				}
				vrrqs[0]->query_start.push_back(vrrqs[0]->seq.size());
				vrrqs[0]->seq+=vrrqs[1]->seq.substr(vrrqs[1]->query_start[i]+shift,
												  vrrqs[1]->block_size[i]-shift);
				vrrqs[0]->target_start.push_back(vrrqs[1]->target_start[i]+shift);
				
				vrrqs[0]->block_size.push_back(vrrqs[1]->block_size[i]-shift);
				vrrqs[0]->size=vrrqs[0]->seq.size();
				vrrqs[0]->is_merged=true;
				// TODO: merge the flag_field
			}
		}
		else { //(vrrqs.size()>3)
			for (i=0;i<vrrqs.size();i++){
				fprintf(stderr,"%s\t%s\t%s\n",vrrqs[i]->name.c_str(),
						vrrqs[i]->target.c_str(),vrrqs[i]->flag_field.c_str());
			}
			fprintf(stderr,
					"WARNING: more than two flags paired read, discarded the read. \n");
			continue;
		}

		// find reference genome, and count mismatches
		for (j=0;j<ref.size();j++){
			if (ref[j]->name==vrrqs[0]->target){
				vrrqs[0]->mismatch=count_mismatches(ref[j],vrrqs[0]);
				break;
			}
		}
		bam_result.push_back(vrrqs[0]);
	}
}

int jeweler::add_queries(vector<Transcript *> &ref, 
				multimap<string,string> &srmap,
				vector<rna_read_query *> &bam_result,
				map<rna_read_key,rna_read_query *>& queries){
	int i,j,k;

	
	// associate each read with its sequence data

	for (i=0;i<bam_result.size();i++){
		
		// double check the read with the aligned segments
		// bam has the exactly same sequence. 
		
		map<string,string>::iterator finder;
		string read_id_in_fasta=bam_result[i]->name+";"+bam_result[i]->flag_field;
		//fprintf(stderr,"%s\n",read_id_in_fasta.c_str());
		if ((finder=srmap.find(read_id_in_fasta))!=srmap.end()){
			bam_result[i]->seq=finder->second;

			for (j=0;j<ref.size();j++){
				if (ref[j]->name==bam_result[i]->target){

					// the number of mismatches does not conform with 
					// the result from bam, which means the seq need to be reversed and 
					// complemented
					int mismatches=count_mismatches(ref[j],bam_result[i]);						

					if (mismatches!=bam_result[i]->mismatch){
						recover_original_read(bam_result[i]->seq);
					}
					mismatches=count_mismatches(ref[j],bam_result[i]);						

					if (mismatches!=bam_result[i]->mismatch)
						printf("%s\t%d\t%d\t%s\n",bam_result[i]->name.c_str(),
							   bam_result[i]->mismatch,mismatches,bam_result[i]->seq.c_str());
					break;
				}
			}
		}
		else {
			fprintf(stderr,"cannot find seq data\n");
			bam_result[i]->is_ignored=true;
			bam_result[i]->is_initialized=false;
			continue;
		}

		// TODO : Does not allow any gap now. 
		int threshold=10;
		int gap_threshold = 15;

		if(bam_result[i]->mismatch>threshold
		   ||bam_result[i]->target_gap_num>gap_threshold
		   || bam_result[i]->matches<(bam_result[i]->size-threshold )
		   ){
			//fprintf(stderr,"%d\t%d\t%d\n",bam_result[i]->mismatch,bam_result[i]->matches,bam_result[i]->size-threshold);
			bam_result[i]->is_ignored=true;
			bam_result[i]->is_initialized=false;
			continue;
		}
		else{
			bam_result[i]->is_ignored=false;
			bam_result[i]->is_initialized=true;
		}

	}
	//fprintf(log_file,"Total %d reads loaded.\n",bam_result.size());
	// refresh bam_result variable with all marged paired read.
	merge_paired_reads(ref,bam_result);
	// insert all bam result into queries
	// may have multiple alignment 
	// choose the best one that can be aligned to either transcript
	//fprintf(log_file,"Total %d reads inserted into reads database.\n",bam_result.size());
	for (i=0;i<bam_result.size();i++){
		rna_read_key rrk=rna_read_key(bam_result[i]->target, bam_result[i]->name);
		map<rna_read_key,rna_read_query *>::iterator iter=queries.find(rrk);

		if (iter==queries.end()){
			queries.insert(make_pair(rrk,bam_result[i]));
		}
		else {
			if (is_better_alignment(bam_result[i],(iter->second))){
				iter->second=bam_result[i];
			}
		}
	}

	return 0;
}

int jeweler::label_mismatches_perbase(std::vector<Transcript*> &ptrans, 
									  std::vector<Transcript*> &mtrans, 
									  std::map<rna_read_key,rna_read_query * > &queries)
{
	std::map<rna_read_key,rna_read_query*>::iterator it;

	for(it=queries.begin(); it !=queries.end(); ++it)
	{

		for(int i = 0; i < ptrans.size(); ++i)
		{
			if (it->second->target==ptrans[i]->name){
				//annotate_mismatch_pos(ptrans[i], it->second);		
				//annotate_mismatch_pos(mtrans[i], it->second);			}
		}
	}
}


int jeweler::load_read_data(TranscriptInfo * ti, 
							vector<Transcript *> &ptrans,
							vector<Transcript *> &mtrans,
							map<rna_read_key , rna_read_query *> &queries
							){
	
	vector<seq_read *> sr;
	vector<rna_read_query *> pqueries,mqueries;
	multimap<string,string> srmap;
	int i;

	//load_fasta_file(ti->read_seq_filename,sr);
	//for (i=0;i<sr.size();i++){
	//srmap.insert(make_pair(sr[i]->name,sr[i]->seq));
	//}
	

	//load_psl_file(ti->paternal_aligned_filename,pqueries);
	//load_psl_file(ti->maternal_aligned_filename,mqueries);
	//add_queries(ptrans,srmap,pqueries,queries);
	//add_queries(mtrans,srmap,mqueries,queries);
	
	// Now, we only need to read one bam file

	return 0;
}

bool jeweler::match_snp(Transcript *t, rna_read_query* rrq){
	if (t->snp_pos.size()==0) 
		return false;
	int i,j;
	int shift;
	int query_allele_pos;
	int no_allele=true;

	for (i=0;i<rrq->target_start.size();i++){
		for (j=0;j<t->snp_pos.size();j++){
			if (t->snp_pos[j]>=rrq->target_start[i]
				&& t->snp_pos[j]<rrq->target_start[i]+rrq->block_size[i]){
				no_allele=false;
				shift=t->snp_pos[j]-rrq->target_start[i];
				query_allele_pos=rrq->query_start[i]+shift;
				if (t->seq[t->snp_pos[j]]!=rrq->seq[query_allele_pos]){
					return false;
				}
			}
		}
	}
	return !no_allele;
}


bool test_match_snp(Transcript * t, Transcript* t2, rna_read_query* rrq){
	if (t->snp_pos.size()==0) 
		return false;
	int i,j;
	int shift;
	int query_allele_pos;
	int no_allele=true;

	for (i=0;i<rrq->target_start.size();i++){
		for (j=0;j<t->snp_pos.size();j++){
			if (t->snp_pos[j]>=rrq->target_start[i]
				&& t->snp_pos[j]<rrq->target_start[i]+rrq->block_size[i]){
				no_allele=false;
				shift=t->snp_pos[j]-rrq->target_start[i];
				query_allele_pos=rrq->query_start[i]+shift;
				if (t->seq[t->snp_pos[j]]!=rrq->seq[query_allele_pos] && 
					t2->seq[t->snp_pos[j]]!=rrq->seq[query_allele_pos]  ){
					printf("%s\t%s\n",rrq->name.c_str(),rrq->seq.c_str());
					printf("%d\t%d\n",rrq->target_start[i],rrq->mismatch);
					printf("%s\t%s\n",t->name.c_str(),t->seq.substr(rrq->target_start[i],rrq->block_size[i]).c_str());
					printf("%s\n",t2->seq.substr(rrq->target_start[i],rrq->block_size[i]).c_str());
				}
			}
		}
	}
	return !no_allele;
}

int jeweler::identify_sources(vector<Transcript *> &source,
							  map<rna_read_key,rna_read_query *> &queries,
							  int source_id){
	int i;

	// because the number of transcripts are no more than five
	// so a full loop should be efficient enough.
	for (i=0;i<source.size();i++){
		map<rna_read_key,rna_read_query *>::iterator j;
		for (j=queries.begin();j!=queries.end();j++){
			if (j->first.transcript_name!=source[i]->name) 
				continue;
			if (match_snp(source[i],j->second)){
				j->second->source_id=source_id;
			}
		}
	}
	return 0;
}
