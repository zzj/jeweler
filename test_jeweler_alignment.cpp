#include <limits.h>
#include <string>
#include "jeweler_alignment.hpp"
#include "gtest/gtest.h"
using std::string;

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
