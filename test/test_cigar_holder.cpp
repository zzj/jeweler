#include "laboratory/cigar_holder.hpp"
#include "main/jeweler_alignment.hpp"
#include "gtest/gtest.h"

TEST(CigarHolderTest, test_get_cigarop) {
    string cigar_string = "10M5N";
    vector<CigarOp> cigar_data;
    get_cigarop(cigar_string, cigar_data);
    ASSERT_EQ(2, cigar_data.size());
    EXPECT_EQ(10, cigar_data[0].Length);
    EXPECT_EQ(5, cigar_data[1].Length);
    EXPECT_EQ('M', cigar_data[0].Type);
    EXPECT_EQ('N', cigar_data[1].Type);
}

TEST(CigarHolderTest, test_get_cigar_string) {
    string cigar_string = "10M5N";
    vector<CigarOp> cigar_data;
    get_cigarop(cigar_string, cigar_data);
    string new_cigar_string = get_cigar_string(cigar_data);
    EXPECT_EQ(cigar_string, new_cigar_string);
}

void test_cigar_trim(unsigned int position, string cigar_string,
                     unsigned int expected_position, string expected_cigar_string){
    JewelerAlignment al;
    al.Position = position;
    get_cigarop(cigar_string, al.CigarData);
    cigar_trim(al);
    string new_cigar_string = get_cigar_string(al);
    EXPECT_EQ(expected_cigar_string, new_cigar_string);
    EXPECT_EQ(expected_position, al.Position);
}


TEST(CigarHolderTest, test_cigar_trim) {
    // test normal one
    test_cigar_trim(1, "10M5N3M", 1, "10M5N3M");
    // test Ns at begining
    test_cigar_trim(1, "5N3M", 6, "3M");
    test_cigar_trim(1, "2N3N3M", 6, "3M");
    // test trailing Ns
    test_cigar_trim(1, "3M2N", 1, "3M");
    test_cigar_trim(1, "3M3N2N", 1, "3M");
    // test both training and begining
    test_cigar_trim(1, "5N3M2N", 6, "3M");
    test_cigar_trim(1, "2N3N3M3N2N", 6, "3M");
    // test Ns in the middle
    test_cigar_trim(1, "3M2N3N5M", 1, "3M5N5M");
    test_cigar_trim(1, "3M2N1N3N5M", 1, "3M6N5M");
    // test complicated cases
    test_cigar_trim(1, "3M2N1N3N5M1N1N3M", 1, "3M6N5M2N3M");
    test_cigar_trim(1, "1N3M2N1N3N5M1N1N3M4N", 2, "3M6N5M2N3M");

}
