#ifndef _TRANSCRIPT_INFO_H_
#define _TRANSCRIPT_INFO_H_
#include <string>
#include <cstdio>

using namespace std;

class transcript_info{
public:

	string gene_id;
	string folder;
	string gtf_filename;
	string paternal_seq_filename;
	string maternal_seq_filename;
	string read_seq_filename;
	
	transcript_info(string gene_id,string folder,string gtf_filename,
					string paternal_seq_filename,
					string maternal_seq_filename,
					string bam_read_filename);

	
};


#endif /* _TRANSCRIPT_INFO_H_ */
