
#include "common.hpp"
using namespace std;

char *trim(char *str)
{
	char *end;

	// Trim leading space
	while(isspace(*str)) str++;

	if(*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace(*end)) end--;

	// Write new null terminator
	*(end+1) = 0;

	return str;
}

FILE * file_open( const char *name, const char * mode) {
     FILE * ret=NULL;
     ret=fopen(name, mode);
     if (ret==NULL)  {
	  fprintf(stderr, "ERROR: can not open file %s  at %s:%d\n",
			     name, __FILE__, __LINE__);
	  exit(0);
     }
     return ret;
}

void open_bam(BamReader &bam_reader, string bam_file){
	if (!bam_reader.Open(bam_file)) {
		fprintf(stderr,"Cannot open bam file %s!\n", bam_file.c_str());
		exit(1);
	}

	if (bam_reader.LocateIndex(BamTools::BamIndex::IndexType::STANDARD) ) {
		fprintf(stdout, "Loaded index for bamfiles ...\n");
	}
	else {
		fprintf(stdout, "Build index for bamfiles ...\n");
		if (!bam_reader.CreateIndex(BamTools::BamIndex::IndexType::STANDARD)) {
			fprintf(stderr,"%s!\n", bam_reader.GetErrorString().c_str());
			exit(1);
		}
	}
}
