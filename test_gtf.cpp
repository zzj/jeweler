#include "gtf.hpp"
#include "transcript.hpp"
#include "gtest/gtest.h"
#include <string>
#include <vector>

using namespace std;

TEST(GTF, test_load_gtf_file) {
	string gtf_filename = "test_data/transcripts.gtf";
	vector<Transcript *> transcripts;
	load_gtf_file(gtf_filename, transcripts);


	EXPECT_EQ(3, transcripts.size());
	EXPECT_EQ(822, transcripts[0]->genome_pos.size());
	EXPECT_EQ(990, transcripts[1]->genome_pos.size());
	EXPECT_EQ(1807, transcripts[2]->genome_pos.size());
	ASSERT_EQ(4, transcripts[0]->exon_start.size());
	ASSERT_EQ(5, transcripts[1]->exon_start.size());
	ASSERT_EQ(2, transcripts[2]->exon_start.size());
	ASSERT_EQ(4, transcripts[0]->exon_end.size());
	ASSERT_EQ(5, transcripts[1]->exon_end.size());
	ASSERT_EQ(2, transcripts[2]->exon_end.size());

    ASSERT_EQ(4, transcripts[0]->num_info_reads_per_exon.size());
    ASSERT_EQ(4, transcripts[0]->num_alleles_per_exon.size());
    ASSERT_EQ(4, transcripts[0]->allele_reads_per_exon.size());
    ASSERT_EQ(4, transcripts[0]->reads_per_exon.size());

    ASSERT_EQ(4766465, transcripts[0]->exon_start[0]);
    ASSERT_EQ(4767606, transcripts[1]->exon_start[1]);
    ASSERT_EQ(4772645, transcripts[2]->exon_start[1]);
    ASSERT_EQ(4766882, transcripts[0]->exon_end[0]);
    ASSERT_EQ(4767729, transcripts[1]->exon_end[1]);
    ASSERT_EQ(4772811, transcripts[2]->exon_end[1]);

	ASSERT_STREQ("CUFF.10", transcripts[2]->gene_id().c_str());
	ASSERT_STREQ("CUFF.9.2", transcripts[1]->transcript_id().c_str());
	ASSERT_STREQ("chr1", transcripts[1]->chr().c_str());
}
