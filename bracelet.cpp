#include "bracelet.hpp"

#include <boost/filesystem.hpp>

using boost::filesystem::path;
using boost::filesystem::create_directory;


Bracelet::Bracelet(JewelerInfo * jeweler_info){
	int id = 0;
	this -> jeweler_info = jeweler_info;
	fprintf(stdout, "bracelet initializing ...\n");

	reads.resize(jeweler_info->gene_id2transcripts.size() );
	reads_index.resize(jeweler_info->gene_id2transcripts.size() );
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
			reads_index[id].insert(temp);
		}
		fclose(foutput_mamf);
		sort(reads[id].begin(),reads[id].end());
		id ++;

	}
}


int Bracelet::intersect(vector<string> &a, vector<string>&b){
	size_t i , j , r;
	i = 0, j = 0;
	r = 0;
	while(i != a.size() && j != b.size()){
		if (a[i]<b[j]) {i++; }
		else if (a[i]>b[j]) {j++; }
		else {i++,j++,r++;}
	}
	return r;
}

int Bracelet::analyze(){
	fprintf(stdout, "bracelet analyzing ...\n");
	for (size_t i=0;i<reads.size(); i++){
		for (size_t j=i+1;j<reads.size();j++){
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

void Bracelet::dump_shared_pileup(FILE * fd,
                                  int original_id,
                                  int target_id){

	map<int, int> coverage;
	char name [1000];
	string fileinput = string(jeweler_info->result_folder + "/" +
                              jeweler_info->gene_id[original_id] +"/" +
                              jeweler_info->gene_id[original_id] +".compatible.reads");
    FILE * finput = fopen(fileinput.c_str(), "r");
    if (finput ==NULL){
        fprintf(stderr, "No compatible reads file supplied : %s\n", fileinput.c_str());
        return ;
    }
    while(fscanf(finput,"%s", name)==1){
        size_t size,t;
        vector<int> temp;

        fscanf(finput, "%zu", &size);
        for(size_t i = 0; i < size; i ++){
            fscanf(finput,"%zu", &t);
            temp.push_back(t);
        }
        if (reads_index[original_id].find(name)==reads_index[original_id].end() ||
            reads_index[target_id].find(name)==reads_index[target_id].end()){
            continue;
        }

        for (size_t i = 0 ; i < temp.size(); i ++){
            if(coverage.find(temp[i])!=coverage.end()){
                coverage[temp[i]]++;
            }
            else{
                coverage[temp[i]]=1;
            }
        }
    }
    fclose(finput);
    for(auto i = coverage.begin(); i != coverage.end(); i ++){
        fprintf(fd,"%d\t%d\n", i->first, i ->second);
    }
}

int Bracelet::dump(FILE * file, string root){
    fprintf(stdout, "bracelet dumping %s...\n", root.c_str());
    for (size_t i=0; i< reads.size() ;i ++){
        if (results[i].size()==0) continue;
        string current_folder = root+"/"+jeweler_info->gene_id[i] + "/";
        create_directory(path(current_folder));

        fprintf(file, "%s\t", jeweler_info->gene_id[i].c_str());
        fprintf(file, "%zu\t", results[i].size());
        for (size_t j = 0; j < results[i].size(); j ++){
            fprintf(file, "%s\t%d\t",
                    jeweler_info->gene_id[related_transcripts[i][j]].c_str(),
                    results[i][j]);
            FILE *fd = fopen((current_folder+ jeweler_info->gene_id[related_transcripts[i][j]]).c_str(), "w+");
            fprintf(stdout, "%s\n", (current_folder+ jeweler_info->gene_id[related_transcripts[i][j]]).c_str());
            dump_shared_pileup(fd, i, related_transcripts[i][j]);

            fclose(fd);
		}
		fprintf(file, "\n");
	}
	return 0;
}
