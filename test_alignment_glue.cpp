#include "alignment_glue.hpp"
#include "gtest/gtest.h"
#include "jeweler_alignment.hpp"
#include "test_jeweler_alignment.hpp"
#include "laboratory/cigar_holder.hpp"

class AlignmentGlueTest : public ::testing::Test {
protected:
	virtual void SetUp() {
    }
	virtual void TearDown() {
    }
    AlignmentGlue ag;
    int skipped_alignment_length;
    vector<CigarOp> cigar_data;

    void test_get_skipped_region(const int pos, const string &cigar_string,
                                 const char rep, const int query_length,
                                 const int test_skipped_length,
                                 const int expected_skipped_length,
                                 const string &expected_string) {
        JewelerAlignment al;
        create_alignment(pos, cigar_string, rep, query_length, al);
        ag.get_skipped_region(&al, test_skipped_length,
                              cigar_data, skipped_alignment_length);
        EXPECT_EQ(expected_skipped_length, skipped_alignment_length);
        EXPECT_EQ(expected_string, get_cigar_string(cigar_data));
    }
};

                             

TEST_F(AlignmentGlueTest, test_get_skipped_region) {
    test_get_skipped_region(10, "10M", '1', 10, 1, 1, "9M");
    test_get_skipped_region(10, "10M", '1', 10, 5, 5, "5M");
    test_get_skipped_region(10, "10M", '1', 10, 10, 10, "");
    test_get_skipped_region(10, "10M", '1', 10, 15, 10, "");
    test_get_skipped_region(10, "10M10N10M", '1', 10, 15, 10, "5N10M");
} 
