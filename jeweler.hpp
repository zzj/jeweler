#ifndef _JEWELER_H_
#define _JEWELER_H_

#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <map>
#include <set>
#include "fasta.hpp"
#include "common.hpp"
#include "gtf.hpp"
#include "rna_read.hpp"
class transcript_info{
public:

	string gene_id;
	string folder;
	string gtf_filename;
	string paternal_seq_filename;
	string maternal_seq_filename;
	string read_seq_filename;
	string paternal_aligned_filename;
	string maternal_aligned_filename;
	
	transcript_info(string gene_id,string folder,string gf,
					string psf,string msf,string tsf,string paf,string maf);
	int load_transcript_info(FILE *);
	
};




class jeweler{
public:
	string info_filename;
	FILE * log_file;
	vector<transcript_info> transcripts_info;
	
	jeweler(int argc, char *argv[]);
	int run();
	int load_info_file();
	int transcript_helper(string read_file,string gtf_file, 
						  vector<transcript> &transcripts);
	int load_transcript_data(transcript_info ti,
							 vector<transcript> &ptrans, 
							 vector<transcript> &mtrans
							 ); // load transcript data 
	int load_read_data(transcript_info ti, map<rna_read_key,rna_read_query>& queires);
	int identify_sources(vector<transcript> source,
						 map<rna_read_key,rna_read_query> &queries,
						 int source_id);

};

#endif /* _JEWELER_H_ */
