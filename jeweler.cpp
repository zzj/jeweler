#include "jeweler.hpp"
#include <cstdio>
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
	int i;
	vector<seq_read> sr;
	load_fasta_file(seq_filename,sr);
	load_gtf_file(gtf_file,transcripts);
	if (transcripts.size()!=sr.size()){
		fprintf(stderr, "ERROR: sequence file (%d sequences) does not match gtf file (%d transcripts) at %s:%d\n",
				(int)sr.size(),(int)transcripts.size(),__FILE__, __LINE__);
		exit(0);

	}
	for (i=0;i<sr.size();i++){
		if (transcripts[i].name != sr[i].name){
			fprintf(stderr, "ERROR: sequence (%s) does not match gtf file (%s)  at %s:%d\n",
					sr[i].name.c_str(),transcripts[i].name.c_str(),__FILE__, __LINE__);
			exit(0);
		}
		transcripts[i].seq=sr[i].seq;
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
		transcript p,m;
		p=ptrans[i];
		m=mtrans[i];
		if (p.seq.size()!=m.seq.size()){

			fprintf(stderr, 
					"ERROR: transcript sequence size does not match at gene %s  at %s:%d\n",
					p.name.c_str(),__FILE__, __LINE__);
			exit(0);
		}

		for (j=0;j<p.seq.size();j++){
			if (p.seq[j]!=m.seq[j]){
				p.snp_pos.push_back(j);
				m.snp_pos.push_back(j);
				p.alleles.push_back(p.seq[j]);
				m.snp_pos.push_back(m.seq[j]);
			}
		}
	}
	return 0;
}
int jeweler::load_read_data(transcript_info ti, 
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
	return 0;

	load_psl_file(ti.paternal_aligned_filename,pqueries);
	load_psl_file(ti.maternal_aligned_filename,mqueries);
	
	// insert all blat result into queries
	// may have multiple alignment 
	// choose the best one that can be aligned to either transcript

	pqueries.insert(pqueries.end(),mqueries.begin(),mqueries.end());
	for (i=0;i<pqueries.size();i++){
		rna_read_key rrk=rna_read_key(pqueries[i].target, pqueries[i].name);
		map<rna_read_key,rna_read_query>::iterator iter=queries.find(rrk);
		if (iter==queries.end())
			queries.insert(make_pair(rrk,pqueries[i]));
		else {
			if (is_better_alignment(pqueries[i],(iter->second))){
				iter->second=pqueries[i];
			}
		}
	}

	// get sequence
	for (map<rna_read_key,rna_read_query>::iterator j=queries.begin();
		 j!=queries.end();
		 j++){
		map<string,string>::iterator iter;
		if ((iter=srmap.find(j->second.name))!=srmap.end()){
			j->second.seq=iter->second;
		}
		else{
			//fprintf(stderr,"Warning, the read %s has no sequence\n",j->second.name.c_str());
			j->second.is_ignored=true;
			j->second.is_initialized=false;
		}
	}
	
	return 0;
}


int jeweler::identify_sources(vector<transcript> source,
							  map<rna_read_key,rna_read_query> &queries,
							  int source_id){
	int i;
	for (i=0;i<source.size();i++){
		
	}
	return 0;
}


int jeweler::run(){
	int i,j;

	vector<transcript> ptrans, mtrans;

	load_info_file();
	for (i=0;i<transcripts_info.size();i++){
		map<rna_read_key,rna_read_query> queries;	
		load_transcript_data(transcripts_info[i],ptrans, mtrans);
		load_read_data(transcripts_info[i],queries);
		identify_sources(ptrans,queries,1);
		identify_sources(ptrans,queries,2);
	}
	return 0;
}

int main(int argc, char * argv[]){
	jeweler j=jeweler(argc, argv);
	j.run();
	return 0;
}
