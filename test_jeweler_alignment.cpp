#include <limits.h>
#include <string>
#include "jeweler_alignment.hpp"
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
