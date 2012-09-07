#include "appraiser.hpp"
using namespace std;

Appraiser::Appraiser(int argc, char* argv[]) {
	int i;
	bool has_info=false;
	log_file=stdout;
	quality_file=stdout;
	test_num=-1;

	for (i=1; i<argc; i++) {

		if (strcmp(argv[i],"-h")==0) {
			fprintf(stdout, "help information for Appraiser:\n");
			fprintf(stdout, " command: appraiser -b bamfile -f fastafile -l logifile -manf mapfile -qf quality file -tn num_of_tests\n");
			exit(0);
		}
		if (strcmp(argv[i],"-b")==0) {
			i++;
			if (i<argc) {
				bamfile=argv[i];
				has_info=true;
			}
		}
		if (strcmp(argv[i],"-f")==0) {
			i++;
			if (i<argc) {
				fastafile=argv[i];
			}
		}
		if (strcmp(argv[i],"-l")==0) {
			i++;
			if (i<argc) {
				log_file=file_open(argv[i],"w+");
			}
		}
		// Multiple Alignment map file
		if (strcmp(argv[i],"-mamf")==0) { 
			i++;
			fprintf(stdout, "here\n");
			if (i<argc) {
				mam_table_file = string(argv[i]) + "_table";
				mam_map_file = string(argv[i]) + "_map";
			}
		}
		//quality file
		if (strcmp(argv[i],"-qf")==0) {
			i++;
			if (i<argc) {
				quality_file=file_open(argv[i],"w+");
			}
		}
		//test number
		if (strcmp(argv[i],"-tn")==0) {
			i++;
			if (i<argc) {
				test_num=atoi(argv[i]);
			}
		}
	}

	if (!has_info) {
		fprintf(stderr, "ERROR: No bam file! You need to specify it by -b\n");
		exit(0);
	}
}

void Appraiser::run() {
	
	BamReader reader;
	JewelerAlignment al;
	Metabam mb;
	SewingMachine sm;
	int i=0;

	if (!reader.Open(bamfile)) {
		fprintf(stderr,"Cannot open bam file!\n");
		exit(0);
	}
	
	mb.initialize(reader);
	sm.initialize(reader, this->mam_map_file);
	while (reader.GetNextAlignment(al)) {
		mb.add_alignment(al);
		sm.add_alignment(al);
		i++;
		if (test_num>0 && i>test_num) break;
	}
	mb.dump_meta_data(log_file);
	mb.dump_mapquality_list(quality_file);
	//sm.output_alignment_connection_map(mam_map_file);
}

int main( int argc, char * argv[]) {
	Appraiser app(argc,argv);
	app.run();
	return 0;
}
