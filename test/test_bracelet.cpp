#include <limits.h>
#include "main/bracelet.hpp"
#include "main/constants.hpp"
#include "gtest/gtest.h"
#include <boost/assign/std/vector.hpp>

using namespace boost::assign; 

using std::string;

class BraceletTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        origin.add_genome_position(100);
        origin.add_genome_position(101);
        origin.add_genome_position(NOT_FOUND);
        origin.add_genome_position(NOT_FOUND);
        target.add_genome_position(NOT_FOUND);
        target.add_genome_position(201);
        target.add_genome_position(202);
        target.add_genome_position(203);
    }

    virtual void TearDown() {
    }

    Jeweler::EarringsData::Read origin, target;
};


TEST_F(BraceletTest, test_add_coverage) {
    map<int, int> coverage;
    add_coverage(origin, coverage);
    ASSERT_EQ(2, coverage.size());
    ASSERT_EQ(1, coverage[100]);
    add_coverage(origin, coverage);
    ASSERT_EQ(2, coverage.size());
    ASSERT_EQ(2, coverage[100]);
    add_coverage(target, coverage);
    ASSERT_EQ(5, coverage.size());
    ASSERT_EQ(2, coverage[100]);
    ASSERT_EQ(1, coverage[201]);
}

TEST_F(BraceletTest, test_add_coverage_details) {
    map<int, map<int, int> > coverage_details;
    add_coverage_details(origin, target, coverage_details);
    ASSERT_EQ(1, coverage_details.size());
    ASSERT_EQ(1, coverage_details[101].size());
    ASSERT_EQ(1, coverage_details[101][201]);
    add_coverage_details(origin, target, coverage_details);
    ASSERT_EQ(1, coverage_details.size());
    ASSERT_EQ(1, coverage_details[101].size());
    ASSERT_EQ(2, coverage_details[101][201]);
    add_coverage_details(origin, origin, coverage_details);
    ASSERT_EQ(2, coverage_details.size());
    ASSERT_EQ(2, coverage_details[101].size());
    ASSERT_EQ(2, coverage_details[101][201]);
    ASSERT_EQ(1, coverage_details[101][101]);

}

TEST_F(BraceletTest, test_get_coverage) {
    map<int, int> coverage;
    add_coverage(origin, coverage);
    EXPECT_FLOAT_EQ(1, get_coverage_rate(origin, coverage));
    EXPECT_FLOAT_EQ(0, get_coverage_rate(target, coverage));
    add_coverage(target, coverage);
    EXPECT_FLOAT_EQ(1, get_coverage_rate(target, coverage));
}
