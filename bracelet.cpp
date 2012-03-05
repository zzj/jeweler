#include "bracelet.hpp"

Bracelet::Bracelet(JewelerInfo * jeweler_info){
	int id = 0;
	this -> jeweler_info = jeweler_info;
	fprintf(stdout, "bracelet initializing ...\n");

	reads.resize(jeweler_info->gene_id2transcripts.size() );
	results.resize(jeweler_info->gene_id2transcripts.size() );
	related_transcripts.resize(jeweler_info->gene_id2transcripts.size() );
	for (auto i = jeweler_info->gene_id2transcripts.begin();
		 i != jeweler_info->gene_id2transcripts.end();
		 i++){
		string gene_id = i->first;
		this->result_folder = jeweler_info->result_folder + "/" +gene_id;
		FILE *foutput_mamf=
			fopen(string(this->result_folder+"/"+gene_id+".mamf.multiple.reads").c_str(),"r");
		if (foutput_mamf == NULL){
			fprintf(stdout,"WARNING: %d\t is not processed by jeweler\n%s\n", id,
					string(this->result_folder+"/"+gene_id+".mamf.multiple.reads").c_str());
			continue;
		}
		if (id % 1 == 0) fprintf(stdout, "%d\n", id);
		char temp[100];
		while(fscanf(foutput_mamf,"%s", temp)==1){
			reads[id].push_back(temp);
		}
		fclose(foutput_mamf);
		sort(reads[id].begin(),reads[id].end());
		id ++;
		
	}
}


int Bracelet::intersect(vector<string> &a, vector<string>&b){
	int i , j , r;
	i = 0, j = 0;
	r = 0;
	while(i!=a.size() && j!=b.size()){
		if (a[i]<b[j]) {i++; }
		else if (a[i]>b[j]) {j++; }
		else {i++,j++,r++;}
	}
	return r;
}

int Bracelet::analyze(){
	fprintf(stdout, "bracelet analyzing ...\n");
	for (int i=0;i<reads.size(); i++){
		for (int j=i+1;j<reads.size();j++){
			int r = intersect(reads[i], reads[j]);
			if (r > 0){
				results[i].push_back(r);
				related_transcripts[i].push_back(j);
				results[j].push_back(r);
				related_transcripts[j].push_back(i);
			}
		}
	}
	return 0;
}

int Bracelet::dump(FILE * file){
	fprintf(stdout, "bracelet dumping ...\n");
	for ( int i=0; i< reads.size() ;i ++){
		fprintf(file, "%s\t", jeweler_info->gene_id[i].c_str());
		fprintf(file, "%d\t", results[i].size());
		for (int j = 0; j < results[i].size(); j ++){
			fprintf(file, "%s\t%d\t",
					jeweler_info->gene_id[related_transcripts[i][j]].c_str(),
					results[i][j]);
		}
		fprintf(file, "\n");
	}
	return 0;
}
