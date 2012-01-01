#include "jeweler.hpp"
#include <cstdio>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

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
	char gene_id[100],folder[100],gtf_filename[100],
		paternal_seq_filename[100],maternal_seq_filename[100],
		bam_read_filename[100];
	int id=0;
	while(fscanf(fd,"%s%s%s%s%s%s",
				 gene_id, folder, gtf_filename,
				 maternal_seq_filename, paternal_seq_filename, bam_read_filename)==6){
		TranscriptInfo * ti= 
			new TranscriptInfo(gene_id,folder,gtf_filename,
								maternal_seq_filename,paternal_seq_filename,
								bam_read_filename);
		transcripts_info.push_back(ti);
		id++;
	}
	fprintf(log_file, "Total %d genes are loaded\n",id);
	return 0;
}




int jeweler::run(){
	int i,j;

	vector<Transcript *> paternal_transcripts, maternal_transcripts;

	load_info_file();
	for (i=0;i<transcripts_info.size();i++){
		if (i%10==0) fprintf(log_file,"%d\n",i);
		Earrings earrings(transcripts_info[i]);
	}
	return 0;
}

int Transcript::tolerate=0;

int main(int argc, char * argv[]){
	jeweler j=jeweler(argc, argv);
	j.run();
	return 0;
}
