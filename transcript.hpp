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
using namespace BamTools;



#include "common.hpp"

using namespace std;
class Transcript{
public:
	Transcript();

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
	template<typename T, typename get_info>
	T get_transcript_aligned_info(BamAlignment *, get_info gi);


	// get aligned sequences from the query.
	string get_query_aligned_seq(BamAlignment *);

    // get transcript locations for a given genome postion
    int get_transcript_location(int genome_location);

	// get covered exons for a given aligment
	int get_transcript_exon(int genome_location);

		
	// get allele char
	char get_allele_char(int transcript_location);
	
	// get the exon given the postion
	int get_allele_exon(int transcript_location);

	// check whether it is allele or not
	bool is_allele(int transcript_location);
	
	// register informative reads by exon
	int register_read(BamAlignment *);

    int output_segments();

	/*****************
	 * Graph utility
	 *****************/
	int add_transcript_to_graph(Graph *);
	
	/******************
     * Landscape plot 
     *****************/
	
	vector<int> num_total_reads; // number of total reads per locus
	int build_landscape_plot();
};


int get_seq_info(Transcript *, int start, int length, string &ret);
int get_exon_info(Transcript *, int start, int length, vector<int>& exons);

template<typename T, typename get_info> 
T Transcript::get_transcript_aligned_info(BamAlignment * al, get_info gi){
	// must perform the compatible test first!
	// justify whether the sequences contains the BamAlignment

	// not sure this is the case, but the coordinates are messed up
	// sometimes. TODO: read the samtools's specification, and make
	// the coordinates right.
	int start_pos=al->Position+1;
	int start;
	T ret;
	
	std::vector< CigarOp > &cigar_data = al->CigarData;
	vector<CigarOp>::const_iterator cigar_iter = cigar_data.begin();
	vector<CigarOp>::const_iterator cigar_end  = cigar_data.end();
	
	for ( ; cigar_iter != cigar_end; ++cigar_iter ) {
		const CigarOp& op = (*cigar_iter);
		
		switch ( op.Type ) {
			
			// for 'M', '=', 'X' - write bases
		case (Constants::BAM_CIGAR_MATCH_CHAR)    :
		case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
		case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
			// the beginning and end of the matched sequence must 
			// be belong to the same sequences. <
			gi(this, start_pos, op.Length, ret);
			start_pos+=op.Length;
			break;
			
		case (Constants::BAM_CIGAR_INS_CHAR)      :			
				// fall through
			break;
		// for 'S' - soft clip, do not write bases
		// but increment placeholder 'k'
		case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :

			break;
			
		// for 'D' - write gap character
		// for 'N' - write N's, skip bases in original query sequence
		// for 'P' - write padding character			
		case (Constants::BAM_CIGAR_DEL_CHAR) :
		case (Constants::BAM_CIGAR_PAD_CHAR) :
		case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
			start_pos+=op.Length;
			break;
			
			// for 'H' - hard clip, do nothing to AlignedBases, move to next op
		case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
			break;
			
			// invalid CIGAR op-code
		default:
			const string message = string("invalid CIGAR operation type: ") + op.Type;
			return false;
		}
	}
	return ret;
	
}

#endif /* _TRANSCRIPT_H_ */
