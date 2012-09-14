#include "gtf.hpp"
#include <boost/algorithm/string.hpp>
#include "fasta.hpp"
#include "common.hpp"
#include "gtf_info.hpp"
#include "constants.hpp"
#include "transcript.hpp"

using namespace std;

void add_new_transcript(const vector<gtf_info>& gtf_list,
                        vector<Transcript *> &maternal_transcripts,
                        vector<Transcript *> &paternal_transcripts,
                        FastaReference * maternal_ref,
                        FastaReference * paternal_ref) {
    Transcript *mt = new Transcript();
    mt->set_origin(TRANSCRIPT_MATERNAL);
    mt->load_gtf(gtf_list);
    mt->load_seq(maternal_ref);

    Transcript *pt = new Transcript();
    pt->set_origin(TRANSCRIPT_PATERNAL);
    pt->load_gtf(gtf_list);
    pt->load_seq(paternal_ref);

    assert(pt->seq().size() == mt->seq().size());
    
    // find SNP position by given both paternal and maternal
    // transcripts sequences
    vector<unsigned int> snp_pos;
    vector<char> maternal_alleles, paternal_alleles;
    for (size_t j = 0; j < pt->seq().size(); j++) {
        if (pt->seq()[j] != mt->seq()[j]) {
            snp_pos.push_back(j);
            maternal_alleles.push_back(mt->seq()[j]);
            paternal_alleles.push_back(pt->seq()[j]);
        }
    }
    
    pt->load_snps(snp_pos, paternal_alleles, EXON_PATERNAL);
    mt->load_snps(snp_pos, maternal_alleles, EXON_MATERNAL);

    maternal_transcripts.push_back(mt);
    paternal_transcripts.push_back(pt);
    return ;
}

// This function is too long, should first get a vector of rows for
// just one transcript 
void load_gtf_file(const string &gtf_filename,
                   vector<Transcript *> &maternal_transcripts,
                   vector<Transcript *> &paternal_transcripts,
                   FastaReference * maternal_ref,
                   FastaReference * paternal_ref) {

	FILE *fd = file_open(gtf_filename.c_str(),"r");
	char * temp = (char *)malloc(MAXLINE+1);
	char * line;

	maternal_transcripts.clear();
	paternal_transcripts.clear();
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
                add_new_transcript(gtf_list,
                                   maternal_transcripts, paternal_transcripts,
                                   maternal_ref, paternal_ref);
                gtf_list.clear();
            }
        }
        gtf_list.push_back(gi);
	}
    add_new_transcript(gtf_list,
                       maternal_transcripts, paternal_transcripts,
                       maternal_ref, paternal_ref);
	free(temp);
	fclose(fd);
}
