#include <limits.h>
#include <string>
#include "main/jeweler_alignment.hpp"
#include "laboratory/cigar_holder.hpp"
#include "test_jeweler_alignment.hpp"
#include "gtest/gtest.h"
using std::string;

void create_alignment(const int position, const string &cigar_string,
                      const char rep, int query_length, JewelerAlignment &al) {
    al.Position = position;
    al.QueryBases = string(query_length, rep);
    al.Qualities = string(query_length, rep);
    get_cigarop(cigar_string, al.CigarData);
}

CreateAlignment::CreateAlignment(): al(new JewelerAlignment) {
}

CreateAlignment& CreateAlignment::set_position(int position) {
    al->Position = position;
    return *this;
}

CreateAlignment& CreateAlignment::set_length(int length) {
    this->length = length;
    return *this;
}

CreateAlignment& CreateAlignment::set_query_bases(char rep) {
    al->QueryBases = string(this->length, rep);
    return *this;
}

CreateAlignment& CreateAlignment::set_qualities(char rep) {
    al->Qualities = string(this->length, rep);
    return *this;
}


CreateAlignment& CreateAlignment::set_cigarop(string cigar_string) {
    get_cigarop(cigar_string, al->CigarData);
    return *this;
}

JewelerAlignment CreateAlignment::get() {
    return *this->al;
}

class JewelerInfoTest : public ::testing::Test {
protected:
	virtual void SetUp() {
		al.QueryBases = string(100, 'G');
	}

	// virtual void TearDown() {}

	JewelerAlignment al;
};

TEST_F(JewelerInfoTest, NormalCase) {
	al.SetIsReverseStrand(false);
	EXPECT_EQ(10, get_read_position(&al, 10));
}

TEST_F(JewelerInfoTest, ReverseCase) {
	al.SetIsReverseStrand(true);
	EXPECT_EQ(89, get_read_position(&al, 10));
}


// TEST AlignmentExpert

class SimpleAlignmentExpert : public AlignmentExpert {
public:
    int len_matched;
    int len_only_read;
    int len_only_genome;
    int len_neither_exist;
    SimpleAlignmentExpert() {
        this->len_matched = 0;
        this->len_only_read = 0;
        this->len_only_genome = 0;
        this->len_neither_exist = 0;
    }
    virtual void study_matched_seq(JewelerAlignment *al, int genome_start,
                                   int alignment_start, int length) {
        this->len_matched += length;
    }

    virtual void study_only_read_seq(JewelerAlignment *al, int genome_start,
                                     int alignment_start, int length) {
        this->len_only_read += length;
    }

    virtual void study_only_genome_seq(JewelerAlignment *al, int genome_start,
                                       int alignment_start, int length) {
        this->len_only_genome += length;
    }

    virtual void study_neither_exist_seq(JewelerAlignment *al, int genome_start,
                                         int alignment_start, int length) {
        this->len_neither_exist += length;
    }
};

TEST(AlignmentExpertTest, SimpleTest) {
    SimpleAlignmentExpert sae;
    JewelerAlignment al= CreateAlignment()
        .set_position(110)
        .set_length(10)
        .set_qualities('T')
        .set_query_bases('G')
        .set_cigarop("1H10M20=30X40N50I60J")
        .get();
    al.investigate(&sae);
    EXPECT_EQ(60, sae.len_matched);
    EXPECT_EQ(100, sae.len_only_genome);
    EXPECT_EQ(50, sae.len_only_read);
    EXPECT_EQ(1, sae.len_neither_exist);
}



class AlignmentGlueTest : public ::testing::Test {
protected:
	virtual void SetUp() {
    }
	virtual void TearDown() {
    }
    int skipped_alignment_length;
    vector<CigarOp> cigar_data;

    void test_get_skipped_region(const int pos, const string &cigar_string,
                                 const char rep, const int query_length,
                                 const int test_skipped_length,
                                 const int expected_skipped_length,
                                 const string &expected_string) {
        JewelerAlignment al;
        create_alignment(pos, cigar_string, rep, query_length, al);
        al.get_skipped_region(test_skipped_length,
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
    test_get_skipped_region(10, "10M10N10M", '1', 10, 25, 15, "5M");
    test_get_skipped_region(10, "5M1D4M10N10M", '1', 10, 15, 9, "5N10M");
    test_get_skipped_region(10, "5M1D4M10N10M", '1', 10, 25, 14, "5M");
    test_get_skipped_region(10, "5M1I4M10N10M", '1', 10, 15, 10, "4N10M");
    test_get_skipped_region(10, "5M1I4M10N10M", '1', 10, 25, 16, "4M");
} 
