#include "gtf.hpp"
#include "transcript.hpp"
#include "gtest/gtest.h"
#include <string>
#include <vector>

using namespace std;

TEST(GTF, test_load_gtf_file) {
	string gtf_filename = "test_data/transcripts.gtf";
    string maternal_strain_ref_file = "test_data/test_transcript_seq_maternal.fa";
    string paternal_strain_ref_file = "test_data/test_transcript_seq_paternal.fa";
	vector<Transcript *> transcripts, paternal_transcripts;

	FastaReference *maternal_fasta = new FastaReference();
	FastaReference *paternal_fasta = new FastaReference();
	maternal_fasta->open(maternal_strain_ref_file);
	paternal_fasta->open(paternal_strain_ref_file);
	load_gtf_file(gtf_filename,
                  transcripts, paternal_transcripts,
                  maternal_fasta, paternal_fasta);

	EXPECT_EQ(3, transcripts.size());
	EXPECT_EQ(320, transcripts[0]->genome_pos.size());
	EXPECT_EQ(240, transcripts[1]->genome_pos.size());
	EXPECT_EQ(160, transcripts[2]->genome_pos.size());
	ASSERT_EQ(4, transcripts[0]->exon_start.size());
	ASSERT_EQ(3, transcripts[1]->exon_start.size());
	ASSERT_EQ(2, transcripts[2]->exon_start.size());
	ASSERT_EQ(4, transcripts[0]->exon_end.size());
	ASSERT_EQ(3, transcripts[1]->exon_end.size());
	ASSERT_EQ(2, transcripts[2]->exon_end.size());

    ASSERT_EQ(4, transcripts[0]->num_info_reads_per_exon.size());
    ASSERT_EQ(4, transcripts[0]->num_alleles_per_exon.size());
    ASSERT_EQ(4, transcripts[0]->allele_reads_per_exon.size());
    ASSERT_EQ(4, transcripts[0]->reads_per_exon.size());

    ASSERT_EQ(161, transcripts[0]->exon_start[0]);
    ASSERT_EQ(481, transcripts[1]->exon_start[1]);
    ASSERT_EQ(641, transcripts[2]->exon_start[1]);
    ASSERT_EQ(240, transcripts[0]->exon_end[0]);
    ASSERT_EQ(560, transcripts[1]->exon_end[1]);
    ASSERT_EQ(720, transcripts[2]->exon_end[1]);

	ASSERT_STREQ("CUFF.10", transcripts[2]->gene_id().c_str());
	ASSERT_STREQ("CUFF.9.2", transcripts[1]->transcript_id().c_str());
	ASSERT_STREQ("transcript1", transcripts[1]->chr().c_str());

    ASSERT_EQ(1, transcripts[0]->snp_pos.size());
    ASSERT_EQ(0, transcripts[0]->snp_pos[0]);
    ASSERT_EQ('C', transcripts[0]->alleles[0]);
    ASSERT_EQ('G', paternal_transcripts[0]->alleles[0]);

}
