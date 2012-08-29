#include <limits.h>
#include "gtest/gtest.h"
#include "../common.hpp"
#include "sewing_machine.hpp"

class SewingMachineTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        al1.Name = "al1";
        al2.Name = "al2";
    }

    virtual void TearDown() {
    }

    void count() {
        sm.count_alignment(al1);
        sm.count_alignment(al2);
        sm.count_alignment(al2);
    }

    void add() {
        sm.add_alignment(al1);
        sm.add_alignment(al2);
        sm.add_alignment(al2);
    }
    SewingMachine sm;
    JewelerAlignment al1;
    JewelerAlignment al2;
};

TEST_F(SewingMachineTest, test_get_full_name) {
    JewelerAlignment al;
    al.Name = "test";
    al.SetIsFirstMate(true);
    EXPECT_STREQ("test\t1", sm.get_full_name(al).c_str());
    al.SetIsFirstMate(false);
    EXPECT_STREQ("test\t2", sm.get_full_name(al).c_str());
}


TEST_F(SewingMachineTest, test_count_alignment) {
    count();
    EXPECT_EQ(2, sm.multiple_alignment_map.size());
    ASSERT_NE(sm.multiple_alignment_map.end(),
              sm.multiple_alignment_map.find(sm.get_full_name(al1)));
    ASSERT_NE(sm.multiple_alignment_map.end(),
              sm.multiple_alignment_map.find(sm.get_full_name(al2)));
    EXPECT_EQ(1, sm.multiple_alignment_map[sm.get_full_name(al1)]);
    EXPECT_EQ(2, sm.multiple_alignment_map[sm.get_full_name(al2)]);
}

TEST_F(SewingMachineTest, test_add_alignment) {
    count();
    add();
    EXPECT_EQ(1, sm.seqs.size());
    ASSERT_EQ(sm.seqs.end(),
              sm.seqs.find(sm.get_full_name(al1)));
    ASSERT_NE(sm.seqs.end(),
              sm.seqs.find(sm.get_full_name(al2)));
    EXPECT_EQ(2, sm.seqs[sm.get_full_name(al2)].size());
}

TEST_F(SewingMachineTest, test_output_and_load) {
    // TODO: build an bam file, output and reload
    // check the consistancy between before and after.
}
