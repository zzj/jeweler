#include "rna_read.hpp"
using namespace std;

rna_read_key::rna_read_key(string tname,string rname){
	transcript_name=tname;
	rna_read_name=rname;
}

bool operator  < (const rna_read_key &a, const rna_read_key &b){
	if (a.transcript_name!=b.transcript_name)
		return a.transcript_name<b.transcript_name;
	else 
		return a.rna_read_name<b.rna_read_name;
}

rna_read_query::rna_read_query(){
	is_initialized=false;
	is_ignored=true;
	is_merged=false;
	source_id=0;
}




bool operator < (const rna_read_query &a, const rna_read_query  &b){
	if (a.name!=b.name) return a.name<b.name;
	else {
		return a.size<b.size;
	}
}

bool is_better_alignment(const rna_read_query &a, const rna_read_query  &b){
	if (a.is_merged!=b.is_merged) return a.is_merged;
	if (a.mismatch!=b.mismatch) return a.mismatch<b.mismatch;
	if (a.target_gap_num!=b.target_gap_num) return a.target_gap_num<b.target_gap_num;
	if (a.matches!=b.matches) return  a.matches<b.matches;
	return true;
}

bool rna_read_query::is_reversed(){
	return this->start_in_query>this->end_in_query;
}
int load_psl_file(string psl_filename, vector<rna_read_query> &queries){


	FILE *fd=file_open(psl_filename.c_str(),"r");
	char * temp=(char *)malloc(MAXLINE+1);
	char * line;

	int i;
	// ingore first five lines
	for (i=0;i<5;i++){
		fgets(temp, MAXLINE, fd);
	}
	while(fgets(temp, MAXLINE, fd)!=NULL){
		rna_read_query q;		
		vector<string> strs;
		vector<string> strs_temp;
		string read_id, flag_field;
		line=trim(temp);
		boost::split(strs,line,boost::is_any_of("\t"));

		q.rep_match=atoi(strs[2].c_str());
		q.unknown_size=atoi(strs[3].c_str());
		q.query_gap_num= atoi(strs[4].c_str());
		q.query_gap_size=atoi(strs[5].c_str());
		q.mismatch=atoi(strs[1].c_str());
		q.matches=atoi(strs[0].c_str());
		q.start_in_query=atoi(strs[11].c_str());
		q.end_in_query=atoi(strs[12].c_str());
		q.target_gap_num=atoi(strs[6].c_str());
		q.target_gap_size=atoi(strs[7].c_str());
		q.target=strs[13];

		// the first part is the id
		// the second part is the flag field, genearted by samtools -x

		boost::split(strs_temp,strs[9],boost::is_any_of(";"));
		if (strs_temp.size()!=2) {
			fprintf(stderr,"the id does not contain flag fileds. please check the first step\n%s\n",strs[9].c_str());
			exit(0);
		}
		read_id=strs_temp[0]; 
		flag_field=strs_temp[1];

		q.name=read_id;
		q.flag_field=flag_field;
		q.size=atoi(strs[10].c_str());

		int num_blocks=atoi(strs[17].c_str());
		boost::split(strs_temp,strs[18],boost::is_any_of(","));
		for (i=0;i<num_blocks;i++){
			q.block_size.push_back(atoi(strs_temp[i].c_str()));
		}
		boost::split(strs_temp,strs[19],boost::is_any_of(","));
		for (i=0;i<num_blocks;i++){
			q.query_start.push_back(atoi(strs_temp[i].c_str()));
		}
		boost::split(strs_temp,strs[20],boost::is_any_of(","));
		for (i=0;i<num_blocks;i++){
			q.target_start.push_back(atoi(strs_temp[i].c_str()));
		}
		queries.push_back(q);

	}

	free(temp);
	fclose(fd);
	return 0;
}

int recover_original_read(string &seq){
	int i;

	for (i=0;i<seq.size();i++){
		if (seq[i]=='A'){
			seq[i]='T';
		}
		else if (seq[i]=='T'){
			seq[i]='A';
		}
		else if (seq[i]=='C'){
			seq[i]='G';
		}
		else if (seq[i]=='G'){
			seq[i]='C';
		}
		else if (seq[i]=='N'){
			seq[i]='N';
		}
		else {
			fprintf(stderr,"Unknown nucleotide!\n");
			fprintf(stderr,"%s\n",seq.c_str());
			exit(0);
		}
	}
	reverse(seq.begin(),seq.end());
}
