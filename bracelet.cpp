#include "bracelet.hpp"

Bracelet::Bracelet(vector<TranscriptInfo *> infos){
	fprintf(stdout, "bracelet initializing ...\n");
	this->info = infos;

	reads.resize(infos.size() );
	results.resize(infos.size() );
	related_transcripts.resize(infos.size() );
	for (int i=0; i< infos.size() ; i++){
		TranscriptInfo * info = infos[i];
		FILE *foutput_mamf=
			fopen(string(info->folder+"/"+info->gene_id+".mamf.multiple.reads").c_str(),"r");
		if (foutput_mamf == NULL){
			fprintf(stdout,"WARNING: %d\t is not processed by jeweler\n%s\n", i,
					string(info->folder+"/"+info->gene_id+".mamf.multiple.reads").c_str());
			continue;
		}
		if (i % 1000 == 0) fprintf(stdout, "%d\n", i);

		char temp[100];
		while(fscanf(foutput_mamf,"%s", temp)==1){
			reads[i].push_back(temp);
		}
		fclose(foutput_mamf);
		sort(reads[i].begin(),reads[i].end());
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
		fprintf(file, "%s\t", info[i]->gene_id.c_str());
		fprintf(file, "%d\t", results[i].size());
		for (int j = 0; j < results[i].size(); j ++){
			fprintf(file, "%s\t%d\t",
					info[related_transcripts[i][j]]->gene_id.c_str(),
					results[i][j]);
		}
		fprintf(file, "\n");
	}
	return 0;
}
