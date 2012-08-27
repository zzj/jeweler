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


	ASSERT_EQ(3, transcripts.size());
	ASSERT_EQ(4, transcripts[0]->exon_start.size());
	ASSERT_EQ(5, transcripts[1]->exon_start.size());
	ASSERT_EQ(2, transcripts[2]->exon_start.size());
	ASSERT_EQ(4, transcripts[0]->exon_end.size());
	ASSERT_EQ(5, transcripts[1]->exon_end.size());
	ASSERT_EQ(2, transcripts[2]->exon_end.size());

	ASSERT_STREQ("CUFF.10", transcripts[2]->gene_id.c_str());
	ASSERT_STREQ("CUFF.9.2", transcripts[1]->transcript_id.c_str());
	ASSERT_STREQ("chr1", transcripts[1]->chr.c_str());

}
