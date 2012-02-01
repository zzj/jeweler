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
#include "constants.hpp"
#include "graph/graph.hpp"
#include "laboratory/cigar_holder.hpp"
#include "alignment_glue.hpp"
using namespace BamTools;



#include "common.hpp"


using namespace std;

class Transcript{

public:
	Transcript();
	
	static int tolerate ;
	bool is_initialized;
	int origin;
	string name;
	string seq;
	string chr;  // storing chromosome name

	// exon info
	// TODO: exon class
	vector<int> exon_start;
	vector<int> exon_end; //end is inclusive 
	vector<int> num_alleles_per_exon; 
	
	// SNP info
	// TODO: snp class
	vector<int> snp_pos;
	vector<char> alleles;
	vector<int> allele_exon; // the number of exon for the given snp

	// the number of informative reads per exon
	vector<int> num_info_reads_per_exon;

	vector<int> genome_pos;  // storing the within chromosome pos of
							 // each base in the transcript
	
	// contains informative reads (with reference alleles)
	// and non informative reads (without reference alleles)
	set<JewelerAlignment *> reads;
	set<JewelerAlignment *> allele_reads;

	vector<set<JewelerAlignment *> > allele_reads_per_exon;
	vector<set<JewelerAlignment *> > reads_per_exon;

	int get_overlapped_alignment(JewelerAlignment *, int &penalty, bool is_to_fix = false);

	// Check the JewelerAlignment is compatible with the transcript
	bool is_compatible(JewelerAlignment *, int tolerate = 0, bool debug = false);

	// Check whether it is a better alignment comparing with other transcript
	int match_alleles(JewelerAlignment *, int &total_alleles,
					  vector<int> &transcript_aligned_locations,
					  vector<int> &read_aligned_locations,
					  vector<int> &alleles, 
					  vector<int> &mismatches,
					  vector<char> &read_mismatch_qualities,
					  vector<char> &mismachars);
	
	// insert the Alignment to the aligned_reads
	int insert_reads(JewelerAlignment *);
	
	// check whether the JewelerAlignment exists in the aligned_reads
	bool is_aligned(JewelerAlignment *);

	// check whether transcripts are the same sequence on the genome
	// or not
	bool is_equal(Transcript *t);


	// get aligned sequences from the transcript.
	template<typename T, typename get_info>
	T get_transcript_aligned_info(JewelerAlignment *, get_info gi);


	// get aligned sequences from the query.
	string get_query_aligned_seq(JewelerAlignment *);

    // get transcript locations for a given genome postion
    int get_transcript_location(int genome_location);

	// get covered exons for a given aligment
	int get_transcript_exon(int genome_location);

	// get next exons
	int get_next_exon(int start_pos, int start_seg , int tolerate);

	// get allele char
	char get_allele_char(int transcript_location);
	
	// get the exon given the postion
	int get_allele_exon(int transcript_location);

	// check whether it is allele or not
	bool is_allele(int transcript_location);
	
	// register reads
	int register_read(JewelerAlignment *);	
	
	// register informative reads by exon
	int register_allele_read(JewelerAlignment *);

    int output_segments();

	/*****************
	 * Graph utility
	 *****************/
	int add_transcript_to_graph(Graph *, vector<Path> &records);
	
	/******************
     * Landscape plot 
     *****************/
	
	vector<int> num_total_reads; // number of total reads per locus
	int build_landscape_plot();
};


int get_seq_info(Transcript *, JewelerAlignment * al, 
				 int genome_start, int alignment_start, int length, 
				 string &ret);

int get_transcript_location_info(Transcript *, JewelerAlignment * al, 
					  int genome_start, int alignment_start, int length, 
					  vector<int>  &ret);
int get_read_location_info(Transcript *, JewelerAlignment * al, 
					  int genome_start, int alignment_start, int length, 
					  vector<int>  &ret);

int get_exon_info(Transcript *,  JewelerAlignment * al, 
				  int genome_start, int alignment_start, int length, 
				  vector<int>& exons);

int insert_mismatch_info(Transcript *,  JewelerAlignment * al, 
						 int start, int length, vector<int> &mismatches);


template<typename T, typename get_info> 
T Transcript::get_transcript_aligned_info(JewelerAlignment * al, get_info gi){
	// must perform the compatible test first!
	// justify whether the sequences contains the JewelerAlignment

	// not sure this is the case, but the coordinates are messed up
	// sometimes. TODO: read the samtools's specification, and make
	// the coordinates right.
	int genome_start=al->Position+1;
	int alignment_start=0;
	T ret;
	
	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - aligned string
		case (Constants::BAM_CIGAR_MATCH_CHAR)    :
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
		case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. <
			gi(this, al,genome_start, alignment_start, op.Length, ret);
			genome_start+=op.Length;
			alignment_start+=op.Length;
			break;
			// none aligned string in reads
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :

			alignment_start+=op.Length;
			break;

			// none aligned string in reference genome
		case 'J'      :						
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
			genome_start+=op.Length;
			break;
			
			// for 'H' - hard clip, do nothing to AlignedBases, move to next op
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			break;
			
			// invalid CIGAR op-code
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			fprintf(stdout, "%s\n", message.c_str());
			exit(0);
			return T();
		}
	}
	return ret;
	
}

#endif /* _TRANSCRIPT_H_ */
