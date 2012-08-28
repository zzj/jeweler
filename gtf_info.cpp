#include "gtf_info.hpp"

using namespace std;

gtf_info::gtf_info(vector<string> & strs) {

	char char_gene_id[100];
	char char_transcript_id[100];

    this->start = atoi(strs[3].c_str());
    this->end = atoi(strs[4].c_str());
    this->chr = strs[0];
    this->type = strs[2];
    
	sscanf(strs[8].c_str(),
           "gene_id \"%[^\"]\"; transcript_id \"%[^\"]\";",
           char_gene_id, char_transcript_id);

    this->gene_id = char_gene_id;
    this->transcript_id = char_transcript_id;
	
}
