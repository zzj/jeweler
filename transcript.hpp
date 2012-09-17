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
#include <Fasta.h>
#include "gtf_info.hpp"
#include "jeweler_alignment.hpp"
#include "common.hpp"

using namespace std;
using namespace BamTools;
class ReadMatcher;
class JewelerAlignment;
class Graph;
class Path;

class Transcript{

private:
    bool is_initialized;

public:
    _PROXY_(Transcript)

    Transcript();
    read_only<string> chr;  // storing chromosome name
    read_only<string> gene_id;
    read_only<string> transcript_id;
    read_only<string> seq;
    static int tolerate ;
    read_only<int> origin;
    read_only<int> start, end;

    // exon info
    // TODO: exon class

    vector<unsigned int> exon_start;
    vector<unsigned int> exon_end; //end is inclusive
    vector<unsigned int> num_alleles_per_exon;

    // SNP info
    // TODO: snp class
    vector<unsigned int> snp_pos;
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

    void load_gtf(const vector<gtf_info> &);

    void load_seq(FastaReference *fr);

    int get_exon_by_genome_pos(unsigned int pos);

    void load_snps(vector<unsigned int> &snp_pos, vector<char> &alleles,
                   int exon_type);

    int get_overlapped_alignment(JewelerAlignment *, int &penalty, bool is_to_fix = false);

    // Check the JewelerAlignment is compatible with the transcript
    bool is_compatible(JewelerAlignment *, int tolerate = 0, bool debug = false);

    // Check whether it is a better alignment comparing with other transcript
    int match_alleles(JewelerAlignment *, int &total_alleles,
                      ReadMatcher *rm);


    // check whether the JewelerAlignment exists in the aligned_reads
    bool is_aligned(JewelerAlignment *);

    // check whether transcripts are the same sequence on the genome
    // or not
    bool is_equal(Transcript *t);

    // get aligned sequences from the query.
    string get_query_aligned_seq(JewelerAlignment *);

    // get transcript locations for a given genome postion
    int get_transcript_location(int genome_location);

    // get covered exons for a given aligment
    int get_transcript_exon(int genome_location);

    // get next exons
    int get_next_exon(int start_pos, size_t start_seg , int tolerate);

    // get allele char
    char get_allele_char(int transcript_location);

    // get the exon given the postion
    int get_allele_exon(int transcript_location);

    // check whether it is allele or not
    bool is_allele(int transcript_location);

    // register reads
    int register_read(JewelerAlignment *);

    // register informative reads by exon
    int register_allele_read(JewelerAlignment *, const ReadMatcher &rm);

    void output_segments();

    /*****************
     * Graph utility
     *****************/
    int add_transcript_to_graph(Graph *, vector<Path> &records);

    void dump_seq(string &result_folder, string &filename);

    /******************
     * Landscape plot
     *****************/

    vector<int> num_total_reads; // number of total reads per locus
    int build_landscape_plot();

    void set_origin(int);
};


#endif /* _TRANSCRIPT_H_ */
