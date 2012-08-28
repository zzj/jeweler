#include <limits.h>
#include "transcript.hpp"
#include "gtf.hpp"
#include "jeweler_info.hpp"
#include "gtest/gtest.h"

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
    string gtf_filename;
    FastaReference *fasta_ref;
	
};

TEST_F(TranscriptTest, test_load_seq) { 
    ASSERT_EQ(320, transcripts[0]->seq.size());
    ASSERT_EQ(320, transcripts[0]->genome_pos.size());
    ASSERT_STREQ("GGTATTTATTTCATTTACATTTCCAATGCTATTCCAAAAGTCCCACACACAGACCCCACTCACTCCCCCACCCACCTAATGCCAACTAGGCCATCTTCTGATACATATGCAGCTATTGAAAGTGCTCTTTGCCTTACAGAAGCTTTGCAATTTTATGAGGCAAGGTTATCCCCACTTTCACCTCTATAAGTTTCAATGTCTATGGTTTTATGTGGAGGTCCTTGATCCACTTGGACTTGAAAAATGCTATTTTTTTTTCTGCTGGATTGTTTTATCTCCTTTGTCAAAGATCAAATGATCATAGGTGTGTGGGTTCATTT", transcripts[0]->seq.c_str());
}

