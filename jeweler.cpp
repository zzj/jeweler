#include "jeweler.hpp"
#include <cstdio>
#include <algorithm>
using namespace std;


transcript_info::transcript_info(string gene_id,string folder,string gf,
								 string psf,string msf,string tsf,string paf,string maf){
	this->gene_id=gene_id;
	this->folder=folder;
	this->gtf_filename=gf;
	this->paternal_seq_filename=psf;
	this->maternal_seq_filename=msf;
	this->read_seq_filename=tsf;
	this->paternal_aligned_filename=paf;
	this->maternal_aligned_filename=maf;
}

jeweler::jeweler(int argc, char * argv[]){
	int i;
	bool has_info=false;
	log_file=stdout;

	for (i=1; i<argc; i++){
		if (strcmp(argv[i],"-i")==0){
			i++;
			if (i<argc){
				info_filename=argv[i];
				has_info=true;
			}
		}
		if (strcmp(argv[i],"-l")==0){
			i++;
			if (i<argc){
				log_file=file_open(argv[i],"w+");
			}
		}
	}

	if (!has_info){
		fprintf(stderr, "ERROR: No infomation file! You need to specify it by -i\n");
		exit(0);
	}
	
}

int jeweler::load_info_file(){
	fprintf(log_file,"Now loading info file ...\n");
	FILE *fd=file_open(info_filename.c_str(),"r");
	char gene_id[100],folder[100],gf[100],psf[100],msf[100],tsf[100],paf[100],maf[100];
	int id=0;
	while(fscanf(fd,"%s%s%s%s%s%s%s%s",gene_id,folder,gf,psf,msf,tsf,paf,maf)==8){
		transcript_info ti=transcript_info(gene_id,folder,gf,psf,msf,tsf,paf,maf);
		transcripts_info.push_back(ti);
		id++;
	}
	fprintf(log_file, "Total %d genes are loaded\n",id);
	return 0;
}

int jeweler::transcript_helper(string seq_filename,string gtf_file, 
							   vector<transcript> &transcripts){
	//assuming gtf file has the same order of transcripts with the seq files
	int i,j;
	vector<seq_read> sr;
	load_fasta_file(seq_filename,sr);
	load_gtf_file(gtf_file,transcripts);
	if (transcripts.size()!=sr.size()){
		fprintf(stderr, "ERROR: sequence file (%d sequences) does not match gtf file (%d transcripts) at %s:%d\n",
				(int)sr.size(),(int)transcripts.size(),__FILE__, __LINE__);
		exit(0);

	}
	for (i=0;i<sr.size();i++){
		for (j=0;j<transcripts.size();j++){
			if (transcripts[j].name != sr[i].name){
				continue;
			}
			transcripts[j].seq=sr[i].seq;
			transcripts[j].noninformative_mismatches.resize(transcripts[j].seq.size(),0);
		}
	}
	return 0;
}
	
int jeweler::load_transcript_data(transcript_info ti,
								  vector<transcript> &ptrans, 
								  vector<transcript> &mtrans
								  ){
	int i,j;

	transcript_helper(ti.paternal_seq_filename,ti.gtf_filename, ptrans );
	transcript_helper(ti.maternal_seq_filename,ti.gtf_filename, mtrans );
	if (ptrans.size()!=mtrans.size() || ptrans.size()==0){
		fprintf(stderr, 
				"ERROR: number of transcripts does not match or no reads at all for gene %s  at %s:%d\n",
				ti.gene_id.c_str(),__FILE__, __LINE__);
		exit(0);
	}

	for (i=0;i<ptrans.size();i++){
		transcript &p=ptrans[i];
		transcript &m=mtrans[i];
		if (p.seq.size()!=m.seq.size()){

			fprintf(stderr, 
					"ERROR: transcript sequence size does not match at gene %s  at %s:%d\n",
					p.name.c_str(),__FILE__, __LINE__);
			exit(0);
		}
		int num=0;
		for (j=0;j<p.seq.size();j++){
			if (p.seq[j]!=m.seq[j]){
				p.snp_pos.push_back(j);
				m.snp_pos.push_back(j);
				p.alleles.push_back(p.seq[j]);
				m.alleles.push_back(m.seq[j]);
				num++;
			}
		}

	}
	return 0;
}

int jeweler::count_mismatches(transcript &t, rna_read_query &rrq){
	if (t.name != rrq.target){
		fprintf(stderr,"transcript names do not match with each other\n");
		exit(0);
	}
	int i,j,mismatches=0;

	for (i=0;i<rrq.block_size.size();i++){
		for (j=0;j<rrq.block_size[i];j++){
			if (rrq.seq[rrq.query_start[i]+j] != t.seq[rrq.target_start[i]+j]){
				mismatches++;
			}
		}
	}
	return mismatches;
}
int jeweler::add_queries(vector<transcript> &ref, 
				multimap<string,string> &srmap,
				vector<rna_read_query> &blat_result,
				map<rna_read_key,rna_read_query>& queries){
	int i,j,k;
	
	// insert all blat result into queries
	// may have multiple alignment 
	// choose the best one that can be aligned to either transcript

	for (i=0;i<blat_result.size();i++){
		rna_read_key rrk=rna_read_key(blat_result[i].target, blat_result[i].name);
		map<rna_read_key,rna_read_query>::iterator iter=queries.find(rrk);


		// TODO : only use the exact one now
		// need to 
		if(blat_result[i].mismatch>10 
		   ||blat_result[i].target_gap_num!=0
		   ||blat_result[i].matches!=blat_result[i].size)
			continue;

		// double check the read with the aligned segments
		// blat has 
		
		map<string,string>::iterator finder;
		if ((finder=srmap.find(blat_result[i].name))!=srmap.end()){
			blat_result[i].seq=finder->second;

			for (j=0;j<ref.size();j++){
				if (ref[j].name==blat_result[i].target){
					int mismatches=count_mismatches(ref[j],blat_result[i]);
					// the number of mismatches does not conform with 
					// the result from blat, means the seq need to be reversed and 
					// complemented


					if (mismatches>blat_result[i].mismatch){
						recover_original_read(blat_result[i].seq);
					}
					mismatches=count_mismatches(ref[j],blat_result[i]);						
					if (mismatches!=blat_result[i].mismatch)
						printf("%s\t%d\t%d\t%s\n",blat_result[i].name.c_str(),blat_result[i].mismatch,mismatches,blat_result[i].seq.c_str());
				}
			}
		}
		else {
			fprintf(stderr,"cannot find seq data\n");
			continue;
		}
		if (iter==queries.end()){
			queries.insert(make_pair(rrk,blat_result[i]));
		}
		else {
			if (is_better_alignment(blat_result[i],(iter->second))){
				iter->second=blat_result[i];
			}
		}
	}


}

int jeweler::label_mismatches_perbase(std::vector<transcript> &ptrans, 
									  std::vector<transcript> &mtrans, 
									  std::map<rna_read_key,rna_read_query> &queries)
{
	std::map<rna_read_key,rna_read_query>::iterator it;
	for(it=queries.begin(); it !=queries.end(); ++it)
	{
		for(int i = 0; i < ptrans.size(); ++i)
		{
			annotate_mismatch_pos(ptrans[i], it->second);		
			annotate_mismatch_pos(mtrans[i], it->second);		
		}
	}
}


int jeweler::load_read_data(transcript_info ti, 
							vector<transcript> &ptrans,
							vector<transcript> &mtrans,
							map<rna_read_key, rna_read_query> &queries
							){
	
	vector<seq_read> sr;
	vector<rna_read_query> pqueries,mqueries;
	multimap<string,string> srmap;
	int i;

	load_fasta_file(ti.read_seq_filename,sr);
	for (i=0;i<sr.size();i++){
		srmap.insert(make_pair(sr[i].name,sr[i].seq));
	}


	load_psl_file(ti.paternal_aligned_filename,pqueries);
	load_psl_file(ti.maternal_aligned_filename,mqueries);
	add_queries(ptrans,srmap,pqueries,queries);
	add_queries(mtrans,srmap,mqueries,queries);
	return 0;
}

bool jeweler::match_snp(transcript t, rna_read_query rrq){
	if (t.snp_pos.size()==0) 
		return false;
	int i,j;
	int shift;
	int query_allele_pos;
	int no_allele=true;

	for (i=0;i<rrq.target_start.size();i++){
		for (j=0;j<t.snp_pos.size();j++){
			if (t.snp_pos[j]>=rrq.target_start[i]
				&& t.snp_pos[j]<rrq.target_start[i]+rrq.block_size[i]){
				no_allele=false;
				shift=t.snp_pos[j]-rrq.target_start[i];
				query_allele_pos=rrq.query_start[i]+shift;
				if (t.seq[t.snp_pos[j]]!=rrq.seq[query_allele_pos]){
					return false;
				}
			}
		}
	}
	return !no_allele;
}

int jeweler::annotate_mismatch_pos(transcript &t, rna_read_query rrq)
{
	int num_mismatches = 0;
	for(int i = 0; i < rrq.query_start.size(); ++i)
	{
		int qstr = rrq.query_start[i];
		int tstr = rrq.target_start[i];
		int len  = rrq.block_size[i];
		for(int j = 0; j < len; ++j)
		{
			if(t.seq[tstr+j]!=rrq.seq[qstr+j]) // find a mismatch or an allele
			{
				num_mismatches++;
				if(find(t.snp_pos.begin(), t.snp_pos.end(), tstr+j)!=t.snp_pos.end())
					continue;
				else
				{
					t.noninformative_mismatches[tstr+j] += 1;
				}
			}
		}
	}
	if (num_mismatches < rrq.mismatch)
	{
		fprintf(stderr, "%s\n", t.seq.c_str());
		fprintf(stderr, "%s\n", rrq.seq.c_str());
		fprintf(stderr, "%d %d \n", num_mismatches, rrq.mismatch);
		fprintf(stderr, "WARNING: inconsistent mismatches at query %s \n", rrq.name.c_str());
	}

	return num_mismatches;

}

bool test_match_snp(transcript t, transcript t2, rna_read_query rrq){
	if (t.snp_pos.size()==0) 
		return false;
	int i,j;
	int shift;
	int query_allele_pos;
	int no_allele=true;

	for (i=0;i<rrq.target_start.size();i++){
		for (j=0;j<t.snp_pos.size();j++){
			if (t.snp_pos[j]>=rrq.target_start[i]
				&& t.snp_pos[j]<rrq.target_start[i]+rrq.block_size[i]){
				no_allele=false;
				shift=t.snp_pos[j]-rrq.target_start[i];
				query_allele_pos=rrq.query_start[i]+shift;
				if (t.seq[t.snp_pos[j]]!=rrq.seq[query_allele_pos] && 
					t2.seq[t.snp_pos[j]]!=rrq.seq[query_allele_pos]  ){
					printf("%s\t%s\n",rrq.name.c_str(),rrq.seq.c_str());
					printf("%d\t%d\n",rrq.target_start[i],rrq.mismatch);
					printf("%s\t%s\n",t.name.c_str(),t.seq.substr(rrq.target_start[i],rrq.block_size[i]).c_str());
					printf("%s\n",t2.seq.substr(rrq.target_start[i],rrq.block_size[i]).c_str());


				}
			}
		}
	}
	return !no_allele;
}

int jeweler::identify_sources(vector<transcript> source,
							  map<rna_read_key,rna_read_query> &queries,
							  int source_id){
	int i;

	// because the number of transcripts are no more than five
	// so a full loop should be efficient enough.
	for (i=0;i<source.size();i++){
		map<rna_read_key,rna_read_query>::iterator j;
		for (j=queries.begin();j!=queries.end();j++){
			if (j->first.transcript_name!=source[i].name) 
				continue;
			if (match_snp(source[i],j->second)){
				j->second.source_id=source_id;
			}
		}
	}
	return 0;
}


int jeweler::generate_landscape(transcript_info ti,
								vector<transcript> &ref,
								map<rna_read_key,rna_read_query> &queries){
	int i,j,k,m;
	string filename=string(ti.folder+ti.gene_id+"landscape.plot.info");
	FILE *fd=file_open(filename.c_str(),"w+");
	for ( i=0;i<ref.size();i++){
		fprintf(stderr,"%d\n",i);
		int size=ref[i].seq.size();
		int num_unknown=0,num_paternal=0,num_maternal=0;
		fprintf(stderr,"%d\n",size);
		vector<int> unknown(size,0);
		vector<int> paternal(size,0);
		vector<int> maternal(size,0);
		vector<int> is_snp(size,0);
		// the last position that a  exon ends
		vector<int> exon_jump(size,0);
		
		for (j=0;j<ref[i].snp_pos.size();j++) {
			is_snp[ref[i].snp_pos[j]]=1;
		}

		int current_pos=0;
		for (j=0;j<ref[i].exon_start.size();j++){
			//exon_end is inclusive
			current_pos=current_pos+ref[i].exon_end[j]-ref[i].exon_start[j]+1;
			exon_jump[current_pos-1]=1;
		}
		map<rna_read_key,rna_read_query>::iterator iter;
		for (iter=queries.begin();iter!=queries.end();iter++){
			if (iter->first.transcript_name==ref[i].name){

				vector<int> *target=NULL;

				if (iter->second.source_id==0){
					num_unknown++;
					target=&unknown;
				}
				if (iter->second.source_id==1){
					num_paternal++;
					target=&paternal;
				}
				if (iter->second.source_id==2){
					num_maternal++;
					target=&maternal;
				}
				if(target==NULL){
					fprintf(stderr,"No source id?!\n");
					exit(0);
				}

				for (k=0;k<iter->second.target_start.size();k++){
					for (m=0;m<iter->second.block_size[k];m++){
						if (iter->second.target_start[k]+m>target->size()){
							fprintf(stderr,"out of range\n");
							exit(0);
						}
						target->at(iter->second.target_start[k]+m)++;
					}
				}
			}
		}		
		for (k=0;k<ref[i].snp_pos.size();k++){
			int target=ref[i].snp_pos[k];
			//printf("%d\t%d\t%d\n",unknown[target],paternal[target],maternal[target]);
		}

		for (k=0;k<ref[i].seq.size();k++){
			int target=k;
			fprintf(fd,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", ref[i].name.c_str(),target,
					unknown[target],paternal[target],maternal[target],
					is_snp[target],exon_jump[target],ref[i].noninformative_mismatches[k]);
		}
	}
	fclose(fd);
	return 0;
}

int jeweler::run(){
	int i,j;

	vector<transcript> ptrans, mtrans;

	load_info_file();
	for (i=0;i<transcripts_info.size();i++){
		if (i%10==0) fprintf(log_file,"%d\n",i);
		map<rna_read_key,rna_read_query> queries;	
		load_transcript_data(transcripts_info[i],ptrans, mtrans);
		load_read_data(transcripts_info[i],ptrans,mtrans,queries);
		identify_sources(ptrans,queries,1);
		identify_sources(mtrans,queries,2);

		label_mismatches_perbase(ptrans, mtrans, queries);



			for (map<rna_read_key,rna_read_query>::iterator j=queries.begin();
				 j!=queries.end();
				 j++){
				//test_match_snp(ptrans[0],mtrans[0],j->second);
				;
			}

		fprintf(log_file,"generating\n");
		generate_landscape(transcripts_info[i],ptrans,queries);
	}
	return 0;
}

int main(int argc, char * argv[]){
	jeweler j=jeweler(argc, argv);
	j.run();
	return 0;
}
