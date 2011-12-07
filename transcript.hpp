#ifndef _TRANSCRIPT_H_
#define _TRANSCRIPT_H_

#include <cstdio>
#include <string> 
#include <cstring> 
#include <cstdlib> 
#include <vector>
#include <algorithm>
#include <set>
#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;



#include "common.hpp"

using namespace std;
class Transcript{
public:
	Transcript();

	bool is_initialized;

	string name;
	string seq;
	string chr;  // storing chromosome name
	// exon info
	// TODO: exon class
	vector<int> exon_start;
	vector<int> exon_end; //end is inclusive 
	
	// SNP info
	// TODO: snp class
	vector<int> snp_pos;
	vector<char> alleles;
	vector<int> allele_exon; // the number of exon for the given snp
	// the number of informative reads per exon
	vector<int> num_info_reads_per_exon;

	vector<int> noninformative_mismatches; // number of non informative mismatches per base
	// number of non informative mismatches equaling to A per base
	vector<int> Anoninformative_mismatches; 
	// number of non informative mismatches equaling to C per base
	vector<int> Cnoninformative_mismatches; 
	// number of non informative mismatches equaling to G per base
	vector<int> Gnoninformative_mismatches; 
	// number of non informative mismatches equaling to T per base
	vector<int> Tnoninformative_mismatches; 
	vector<int> genome_pos;  // storing the within chromosome pos of
							 // each base in the transcript
	
	// contains informative reads (with reference alleles)
	// and non informative reads (without reference alleles)
	set<BamAlignment *> aligned_reads;
	set<BamAlignment *> allele_aligned_reads;

	// Check the BamAlignment is compatible with the transcript
	bool is_compatible(BamAlignment *);

	// Check whether it is a better alignment comparing with other transcript
	int match_alleles(BamAlignment *, int &total_alleles, int &match_alleles,
					  vector<int> &alleles);
	
	// insert the Alignment to the aligned_reads
	int insert_aligned_reads(BamAlignment *);

	// check whether the BamAlignment exists in the aligned_reads
	bool is_aligned(BamAlignment *);

	// get aligned sequences from the transcript.
	string get_transcript_aligned_seq(BamAlignment *);

	// get aligned sequences from the query.
	string get_query_aligned_seq(BamAlignment *);

    // get transcript locations for a given genome postion
    int get_transcript_location(int genome_location);

	// get allele char
	char get_allele_char(int transcript_location);
	
	// get the exon given the postion
	int get_allele_exon(int transcript_location);

	// check whether it is allele or not
	bool is_allele(int transcript_location);
	
	// register informative reads by exon
	int register_read(BamAlignment *);

    int output_segments();
	
};


#endif /* _TRANSCRIPT_H_ */
