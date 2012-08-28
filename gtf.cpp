#include "gtf.hpp"
#include <boost/algorithm/string.hpp>
using namespace std;

// This function is too long, should first get a vector of rows for
// just one transcript 
void load_gtf_file(const string &gtf_filename, vector<Transcript *> &transcripts) {

	FILE *fd = file_open(gtf_filename.c_str(),"r");
	char * temp = (char *)malloc(MAXLINE+1);
	char * line;
	Transcript* t;

	transcripts.clear();
    vector<gtf_info> gtf_list;

	while(fgets(temp, MAXLINE, fd)!=NULL) {
		vector<string> strs;
		line = trim(temp);
		boost::split(strs, line, boost::is_any_of("\t"));
        gtf_info gi(strs);
        if (strs.size() != 9) {
            fprintf(stderr, "ERROR: GTF file %s corrupted \n",
                    gtf_filename.c_str());
            exit(0);
        }

        if (gi.type == "transcript") {
            if (gtf_list.size() != 0) {
                t = new Transcript();
                t->load_gtf(gtf_list);
                transcripts.push_back(t);
                gtf_list.clear();
            }
        }
        gtf_list.push_back(gi);
	}
    t = new Transcript();
    t->load_gtf(gtf_list);
	transcripts.push_back(t);
	free(temp);
	fclose(fd);
}
