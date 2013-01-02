#include "gtest/gtest.h"
#include "main/jeweler_info.hpp"

using namespace std;

TEST(JelewerInfoTest, test_check_args) {
	char * argv[3];
	argv[0] = "./jeweler";
	argv[1] = "-bam_file";
	argv[2] = "bam.file";
	string a;
	JewelerInfo ji;
	ASSERT_EQ(2, ji.check_args(1, argv, "-bam_file", a));
	ASSERT_STREQ("bam.file", a.c_str());
}
