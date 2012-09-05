#include "jeweler.hpp"
#include <cstdio>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "proto/jeweler.pb.h"

using namespace std;

jeweler::jeweler(int argc, char * argv[]) {
	
	int i;
	bool has_info = false;
	is_earrings = false;
	is_bracelet = false;
	is_mismatch_analyzer = false;
	is_prepare = false;
	test_case = -1;
	sm = NULL;
	log_file = stdout;
	mamf_filename = "none";
	jeweler_info = new JewelerInfo(argc, argv);
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-h") == 0) {
			fprintf(stdout, "help information for jeweler:\n");
			fprintf(stdout, "command: jeweler -i infofile [-l logfile] [-mamf multiple_alignment_file] -earrnings -bracelet -mismatch_analyzer\n");
			exit(0);
		}
		if (strcmp(argv[i], "-i") == 0) {
			i++;
			if (i < argc) {
				info_filename = argv[i];
				has_info = true;
			}
		}
		if (strcmp(argv[i], "-mamf") == 0) {
			i++;
			if (i < argc) {
				mamf_filename = argv[i];
			}
		}
		if (strcmp(argv[i], "-l") == 0) {
			i++;
			if (i < argc) {

				log_file = file_open(argv[i],"w+");
			}
		}
		if (strcmp(argv[i], "-bracelet") == 0) {
			is_bracelet = true;
			i ++;
			if (i < argc) {
				bracelet_filename = argv[i];
			}
		}
		if (strcmp(argv[i], "-earrings") == 0) {
			is_earrings = true;
		}
		if (strcmp(argv[i], "-prepare") == 0) {
			is_prepare = true;
		}
		if (strcmp(argv[i], "-mismatch_analyzer") == 0) {
			is_mismatch_analyzer = true;
			i ++;
			if (i < argc) {
				mismatch_filename = argv[i];
			}
		}
		if (strcmp(argv[i], "-testcase") == 0) {
			i++;
			if (i < argc) {
				sscanf(argv[i],"%d", &test_case);
			}
		}
	}

	if (!has_info) {
		fprintf(stderr, "ERROR: No infomation file! You need to specify it by -i\n");
		exit(0);
	}
	
}

void jeweler::load_mamf_file() {
	if (mamf_filename != "none") {
		fprintf(stdout,"loading mamf files ...\n");
		sm =new SewingMachine();
		sm->load_multiple_alignments_set(mamf_filename);
	}
}

int jeweler::run() {

	if (is_earrings) {

		load_mamf_file();
		int id = -1;
		for (auto i = jeweler_info->gene_id2transcripts.begin();
			 i != jeweler_info->gene_id2transcripts.end();
			 i++) {
			if (i->second[0]->chr =="chrM") continue;
			create_directory(path(jeweler_info->result_folder + "/" + i->first));
			id ++;
			if (test_case > 0  && id != test_case)  continue;
			else {
				if (test_case > 0) {

				}
			}
			if (id%10 == 0) fprintf(log_file,"%d\n",id);
			printf("%s\n", i->first.c_str());
			Earrings earrings(jeweler_info, 
							  i->first,
							  sm,
							  is_prepare);
		}
	}
	if (!is_prepare) {
		if (is_bracelet) {
			fprintf(stdout, "Bracelet analyzing ...\n");
			Bracelet bracelet(jeweler_info);
			bracelet.analyze();
			FILE * output = fopen((bracelet_filename+"/result.bracelet").c_str(), "w+");
			if (output == NULL) {
				fprintf(stdout, "cannot open file %s\n", (bracelet_filename+"/result.bracelet").c_str());
				exit(0);
			}
			bracelet.dump(output, bracelet_filename);
			fclose(output);
		}
	
		if (is_mismatch_analyzer) {
			fprintf(stdout, "Mismatch analyzing ...\n");
			TranscriptMismatcherAnalyzer tma(mismatch_filename, jeweler_info);
			tma.analyze();
		}
	}
	return 0;
}

int main(int argc, char * argv[]) {
    GOOGLE_PROTOBUF_VERIFY_VERSION;

	jeweler j=jeweler(argc, argv);
	j.run();
	return 0;
}
