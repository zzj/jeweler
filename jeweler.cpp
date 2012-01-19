#include "jeweler.hpp"
#include <cstdio>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

jeweler::jeweler(int argc, char * argv[]){
	
	int i;
	bool has_info = false;
	is_earrings = false;
	is_bracelet = false;
	test_case = -1;
	sm = NULL;
	log_file = stdout;
	mamf_filename = "none";
	for (i = 1; i < argc; i++){
		if (strcmp( argv[i], "-h") == 0){
			fprintf( stdout, "help information for jeweler:\n");
			fprintf( stdout, "command: jeweler -i infofile [-l logfile] [-mamf multiple_alignment_file] \n");
			exit(0);
		}
		if (strcmp( argv[i], "-i") == 0){
			i++;
			if (i < argc){
				info_filename = argv[i];
				has_info = true;
			}
		}
		if (strcmp( argv[i], "-mamf") == 0){
			i++;
			if (i < argc){
				mamf_filename = argv[i];
			}
		}
		if (strcmp( argv[i], "-l") == 0){
			i++;
			if (i < argc){
				log_file = file_open(argv[i],"w+");
			}
		}
		if (strcmp( argv[i], "-bracelet") == 0){
			is_bracelet = true;
		}
		if (strcmp( argv[i], "-earrings") == 0){
			is_earrings = true;
		}
		if (strcmp( argv[i], "-testcase") == 0){
			i++;
			if (i < argc){
				sscanf(argv[i],"%d", &test_case );
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
	char gene_id[300],folder[300],gtf_filename[300],
		paternal_seq_filename[300],maternal_seq_filename[300],
		bam_read_filename[300];
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
		//fprintf(stdout, "%d\n",id);
	}
	fprintf(log_file, "Total %d genes are loaded\n",id);
	fclose(fd);
	return 0;
}


int jeweler::load_mamf_file(){
	if (mamf_filename != "none") {
		FILE * fd = fopen(mamf_filename.c_str(), "r");
		sm =new SewingMachine();
		sm->load_multiple_alignments_set(fd);
		fclose(fd);
	}
}

int jeweler::run(){
	int i,j;

	vector<Transcript *> paternal_transcripts, maternal_transcripts;
	load_info_file();
	if (is_earrings) {

		load_mamf_file();
		if (test_case > 0 ) {
			Earrings earrings(transcripts_info[test_case], sm);
		}
		else {
			for (i = 0; i < transcripts_info.size(); i++){
				if (i%10 == 0) fprintf(log_file,"%d\n",i);
				Earrings earrings(transcripts_info[i], sm);
			}
		}
	}
	if (is_bracelet){
		FILE * output = fopen((info_filename+".bracelet").c_str(), "w+");
		if (output == NULL){
			fprintf(stdout, "cannot open file %s\n", (info_filename+".bracelet").c_str());
			exit(0);
		}
		Bracelet bracelet(transcripts_info);
		bracelet.analyze();
		bracelet.dump(output);
		fclose(output);
	}
	return 0;
}

int Transcript::tolerate=0;

int main(int argc, char * argv[]){
	jeweler j=jeweler(argc, argv);
	j.run();
	return 0;
}
