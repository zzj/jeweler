#include "jeweler.hpp"
#include <cstdio>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

int JewelerInfo::check_args(int &i, char *argv[], const char * name, string &a){
	if (strcmp(argv[i], name) == 0){
		i++;
		a = argv[i];
	}
	return 0;
}


int JewelerInfo::build_gene_id2transcripts(){
	fprintf(stdout, "Bulding gene_id to transcripts index ...\n");
	int i;
	gene_id2transcripts.clear();
	for ( i = 0; i < transcripts.size(); i++){
		gene_id2transcripts[transcripts[i]->gene_id].push_back(transcripts[i]);
	}
	gene_id.clear();
	for ( auto j = gene_id2transcripts.begin();
		  j != gene_id2transcripts.end();
		  j ++){
		gene_id.push_back(j->first);
	}
	return 0;
}

int JewelerInfo::get_refID(string chr){
	for ( int i = 0; i < references.size(); i++){
		if (references[i].RefName == chr){
			return i;
		}
	}
	return NOT_FOUND;
}

JewelerInfo::JewelerInfo(int argc, char *argv []){
	int i;
	for ( i = 1; i < argc; i ++){
		check_args(i, argv, "-maternal_strain_ref_file", maternal_strain_ref_file);
		check_args(i, argv, "-paternal_strain_ref_file", paternal_strain_ref_file);
		check_args(i, argv, "-maternal_strain_id", maternal_strain_id);
		check_args(i, argv, "-paternal_strain_id", paternal_strain_id);
		check_args(i, argv, "-bam_file", bam_file);
		check_args(i, argv, "-alias", alias);
		check_args(i, argv, "-result_file", result_file);
		check_args(i, argv, "-result_folder", result_folder);
		check_args(i, argv, "-gtf_input_file", gtf_input_file);
	}
	result_folder += string("/") + alias + "/";
	create_directory(path(result_folder));
	load_gtf_file(gtf_input_file, transcripts);
	build_gene_id2transcripts();
	
	fprintf(stdout, "There are totally %d transcripts loaded from %s\n", 
			transcripts.size(), gtf_input_file.c_str());
	maternal_fasta = new FastaReference();
	paternal_fasta = new FastaReference();
	maternal_fasta->open(maternal_strain_ref_file);
	paternal_fasta->open(paternal_strain_ref_file);

	if (!bam_reader.Open(bam_file)){
		fprintf(stderr,"Cannot open bam file %s!\n", bam_file.c_str());
		exit(0);
	}

	if (bam_reader.LocateIndex(BamTools::BamIndex::IndexType::STANDARD) ) {
		fprintf(stdout, "Loaded index for bamfiles ...\n");
	}
	else {
		fprintf(stdout, "Build index for bamfiles ...\n");
		if (!bam_reader.CreateIndex(BamTools::BamIndex::IndexType::STANDARD)){
			fprintf(stderr,"%s!\n", bam_reader.GetErrorString().c_str());
			exit(0);
		}
	}
	references = bam_reader.GetReferenceData();
}

JewelerInfo::~JewelerInfo(){

}

jeweler::jeweler(int argc, char * argv[]) {
	
	int i;
	bool has_info = false;
	is_earrings = false;
	is_bracelet = false;
	is_mismatch_analyzer = false;
	test_case = -1;
	sm = NULL;
	log_file = stdout;
	mamf_filename = "none";
	jeweler_info = new JewelerInfo(argc, argv);
	for (i = 1; i < argc; i++){
		if (strcmp( argv[i], "-h") == 0){
			fprintf( stdout, "help information for jeweler:\n");
			fprintf( stdout, "command: jeweler -i infofile [-l logfile] [-mamf multiple_alignment_file] -earrnings -bracelet -mismatch_analyzer\n");
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
			i ++;
			if (i < argc) {
				bracelet_filename = argv[i];
			}
		}
		if (strcmp( argv[i], "-earrings") == 0){
			is_earrings = true;
		}
		if (strcmp( argv[i], "-mismatch_analyzer") == 0){
			is_mismatch_analyzer = true;
			i ++;
			if (i < argc){
				mismatch_filename = argv[i];
			}
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

	if (is_earrings) {

		load_mamf_file();
		int id = 0;
		for (auto i = jeweler_info->gene_id2transcripts.begin();
			 i != jeweler_info->gene_id2transcripts.end();
			 i++){
			if (i->second[0]->chr =="chrM") continue;
			create_directory(path(jeweler_info->result_folder + "/" + i->first));
			if (test_case > 0  && id != test_case)  continue;
			if (id%10 == 0) fprintf(log_file,"%d\n",id);
			id ++;
			Earrings earrings(jeweler_info, 
							  i->first,
							  sm);
		}
	}
	
	if (is_bracelet){
		fprintf(stdout, "Bracelet analyzing ...\n");
		Bracelet bracelet(jeweler_info);
		bracelet.analyze();
		FILE * output = fopen((bracelet_filename+".bracelet").c_str(), "w+");
		if (output == NULL){
			fprintf(stdout, "cannot open file %s\n", (bracelet_filename+".bracelet").c_str());
			exit(0);
		}
		bracelet.dump(output);
		fclose(output);
	}
	
	if (is_mismatch_analyzer){
		fprintf(stdout, "Mismatch analyzing ...\n");
		TranscriptMismatcherAnalyzer tma(mismatch_filename, jeweler_info);
		tma.analyze();
	}
	return 0;
}

int Transcript::tolerate=0;

int main(int argc, char * argv[]){
	jeweler j=jeweler(argc, argv);
	j.run();
	return 0;
}

