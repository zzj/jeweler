#ifndef _EARRINGS_H_
#define _EARRINGS_H_



#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <Fasta.h>
#include "jeweler.hpp"
#include "transcript_info.hpp"
using namespace BamTools;
using namespace std;

// The earrings always present as a pair, though sometimes you may
// only find one. 

// This class will investigate the possible allele specific pile up graph.


class Earrings{
public :
	Earrings(transcript_info *);


	transcript_info * info;
	vector<transcript *> maternal_transcripts;
	vector<transcript *> paternal_transcripts;

	
	// create paternal and maternal transcripts databases
	// annotate SNPs for paternal  and maternal transcripts' sequence
	int load_transcript_data(transcript_info * ti); // load transcript data 

	// load cufflinks' gtf file 
	int transcript_helper(string read_file,string gtf_file, 
						  vector<transcript *> &transcripts);



	
};


#endif /* _EARRINGS_H_ */
