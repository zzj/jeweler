#include <limits.h>
#include "transcript.hpp"
#include "gtf.hpp"
#include "jeweler_info.hpp"
#include "gtest/gtest.h"
#include "test.hpp"
#include "laboratory/cigar_holder.hpp"
#include <boost/assign/std/vector.hpp>

using namespace boost::assign; 

using std::string;

class TranscriptTest : public ::testing::Test {
protected:
	virtual void SetUp() {
        this->gtf_filename = "test_data/test_transcript.gtf";
        string fasta_name = "test_data/test_genome.fa";
        
        this->fasta_ref = new FastaReference();
        this->fasta_ref->open(fasta_name);
        load_gtf_file(gtf_filename, this->transcripts);
        this->transcripts[0]->load_seq(this->fasta_ref);
    }
    
	virtual void TearDown() {
        for (size_t i = 0; i < transcripts.size(); i ++) {
            delete transcripts[i];
        }
    }
    vector<Transcript *> transcripts;
    vector<Transcript *> paternal_transcripts;
    vector<Transcript *> maternal_transcripts;
    string gtf_filename;
    FastaReference *fasta_ref;
};

TEST_F(TranscriptTest, test_load_seq) { 
    ASSERT_EQ(320, transcripts[0]->seq.size());
    ASSERT_EQ(320, transcripts[0]->genome_pos.size());
    ASSERT_STREQ("GGTATTTATTTCATTTACATTTCCAATGCTATTCCAAAAGTCCCACACACAGACCCCACTCACTCCCCCACCCACCTAATGCCAACTAGGCCATCTTCTGATACATATGCAGCTATTGAAAGTGCTCTTTGCCTTACAGAAGCTTTGCAATTTTATGAGGCAAGGTTATCCCCACTTTCACCTCTATAAGTTTCAATGTCTATGGTTTTATGTGGAGGTCCTTGATCCACTTGGACTTGAAAAATGCTATTTTTTTTTCTGCTGGATTGTTTTATCTCCTTTGTCAAAGATCAAATGATCATAGGTGTGTGGGTTCATTT", transcripts[0]->seq.c_str());
}

TEST_F(TranscriptTest, test_load_gtf) {
    // The function is tested in test_gtf.cpp
}

TEST_F(TranscriptTest, test_load_snps) {
    vector<unsigned int> snp_pos;
    snp_pos += 0, 1, 81, 161, 241;
    vector<int> allele_exon;
    allele_exon += 0, 0, 1, 2, 3;
    vector<int> num_alleles_per_exon;
    num_alleles_per_exon += 2, 1, 1, 1;
    vector<char> alleles;
    alleles += 'A', 'T', 'C', 'G';

    transcripts[0]->load_snps(snp_pos, alleles, EXON_PATERNAL);
    EXPECT_EQ(snp_pos.size(), transcripts[0]->snp_pos.size());
    EXPECT_EQ(alleles.size(), transcripts[0]->alleles.size());
    EXPECT_ITERABLE_EQ(vector<unsigned int>, snp_pos, transcripts[0]->snp_pos);
    EXPECT_ITERABLE_EQ(vector<int>, allele_exon, transcripts[0]->allele_exon);
    EXPECT_ITERABLE_EQ(vector<int>, num_alleles_per_exon,
                       transcripts[0]->num_alleles_per_exon);
    EXPECT_ITERABLE_EQ(vector<char>, alleles, transcripts[0]->alleles);
    EXPECT_EQ(EXON_PATERNAL, transcripts[0]->origin);
    
}

TEST_F(TranscriptTest, test_get_exon_by_genome_pos) {
    unsigned int pos;
    pos = transcripts[0]->genome_pos[0];
    EXPECT_EQ(0, transcripts[0]->get_exon_by_genome_pos(pos));
    pos = transcripts[0]->genome_pos[319];
    EXPECT_EQ(3, transcripts[0]->get_exon_by_genome_pos(pos));
    pos = transcripts[0]->genome_pos[319] + 1;
    EXPECT_EQ(-1, transcripts[0]->get_exon_by_genome_pos(pos));
}

TEST_F(TranscriptTest, test_compatible) {
    BamReader bam_reader;
    vector<JewelerAlignment *> compatible_reads;
    open_bam(bam_reader, "test_data/tophat_out/accepted_hits.bam");
    JewelerAlignment *al=new JewelerAlignment();
    while(bam_reader.GetNextAlignment(*al)) {
        compatible_reads.push_back(al);
        al=new JewelerAlignment();
    }
    return ;
    for (size_t i = 0; i < compatible_reads.size(); i++) {
        int32_t ed = 0; // edit distance
        EXPECT_EQ(true, transcripts[0]->is_compatible(compatible_reads[i]));
        fprintf(stderr, "%zu\n", i);
        fprintf(stderr, "%d\n", compatible_reads[i]->Position);
        fprintf(stderr, "%s\n", get_cigar_string(*compatible_reads[i]).c_str());
        fprintf(stderr, "%s\n", compatible_reads[i]->QueryBases.c_str());
    }

    for (size_t i = 0; i < compatible_reads.size(); i++) {
        delete compatible_reads[i];
    }
 }
