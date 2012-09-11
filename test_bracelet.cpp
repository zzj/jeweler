#include <limits.h>
#include "bracelet.hpp"
#include "gtest/gtest.h"
#include <boost/assign/std/vector.hpp>

using namespace boost::assign; 

using std::string;

class BraceletTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        origin.add_read_position(0);
        origin.add_genome_position(0);

        origin.add_read_position(1);
        target.add_read_position(1);
        origin.add_genome_position(1);
        target.add_genome_position(2);

        target.add_read_position(4);
        target.add_genome_position(5);
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
    ASSERT_EQ(2, coverage[1]);
}

TEST_F(BraceletTest, test_add_coverage_details) {
    map<int, map<int, int> > coverage_details;
    add_coverage_details(origin, target, coverage_details);
    ASSERT_EQ(1, coverage_details.size());
    ASSERT_EQ(1, coverage_details[1].size());
    ASSERT_EQ(1, coverage_details[1][2]);
    add_coverage_details(origin, target, coverage_details);
    ASSERT_EQ(1, coverage_details.size());
    ASSERT_EQ(1, coverage_details[1].size());
    ASSERT_EQ(2, coverage_details[1][2]);
    add_coverage_details(origin, origin, coverage_details);
    ASSERT_EQ(2, coverage_details.size());
    ASSERT_EQ(2, coverage_details[1].size());
    ASSERT_EQ(2, coverage_details[1][2]);
    ASSERT_EQ(1, coverage_details[1][1]);
}

TEST_F(BraceletTest, test_get_coverage) {
    map<int, int> coverage;
    add_coverage(origin, coverage);
    ASSERT_FLOAT_EQ(1, get_coverage_rate(origin, coverage));
    ASSERT_FLOAT_EQ(0, get_coverage_rate(target, coverage));
    origin.add_genome_position(3);
    origin.add_genome_position(4);
    ASSERT_FLOAT_EQ(0.5, get_coverage_rate(origin, coverage));
}
