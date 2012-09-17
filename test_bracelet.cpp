#include <limits.h>
#include "bracelet.hpp"
#include "constants.hpp"
#include "gtest/gtest.h"
#include <boost/assign/std/vector.hpp>

using namespace boost::assign; 

using std::string;

class BraceletTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        origin.add_genome_position(0);
        origin.add_genome_position(1);
        origin.add_genome_position(NOT_FOUND);
        origin.add_genome_position(NOT_FOUND);

        target.add_genome_position(NOT_FOUND);
        target.add_genome_position(1);
        target.add_genome_position(2);
        target.add_genome_position(3);
    }

    virtual void TearDown() {
    }

    Jeweler::EarringsData::Read origin, target;
};


TEST_F(BraceletTest, test_add_coverage) {
    map<int, int> coverage;
    add_coverage(origin, coverage);
    ASSERT_EQ(2, coverage.size());
    ASSERT_EQ(1, coverage[1]);
    add_coverage(origin, coverage);
    ASSERT_EQ(2, coverage.size());
    ASSERT_EQ(2, coverage[1]);
    add_coverage(target, coverage);
    ASSERT_EQ(4, coverage.size());
    ASSERT_EQ(3, coverage[1]);
}

TEST_F(BraceletTest, test_add_coverage_details) {
    map<int, map<int, int> > coverage_details;
    add_coverage_details(origin, target, coverage_details);
    ASSERT_EQ(1, coverage_details.size());
    ASSERT_EQ(1, coverage_details[1].size());
    ASSERT_EQ(1, coverage_details[1][1]);
    add_coverage_details(origin, target, coverage_details);
    ASSERT_EQ(1, coverage_details.size());
    ASSERT_EQ(1, coverage_details[1].size());
    ASSERT_EQ(2, coverage_details[1][1]);
    add_coverage_details(origin, origin, coverage_details);
    ASSERT_EQ(2, coverage_details.size());
    ASSERT_EQ(1, coverage_details[1].size());
    ASSERT_EQ(3, coverage_details[1][1]);

}

TEST_F(BraceletTest, test_get_coverage) {
    map<int, int> coverage;
    add_coverage(origin, coverage);
    EXPECT_FLOAT_EQ(1, get_coverage_rate(origin, coverage));
    EXPECT_FLOAT_EQ(1.0/3, get_coverage_rate(target, coverage));
}
