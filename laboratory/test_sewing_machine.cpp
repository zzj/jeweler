#include <limits.h>
#include "gtest/gtest.h"
#include "../common.hpp"
#include "sewing_machine.hpp"
#include <api/BamReader.h>
#include <boost/filesystem.hpp>

using boost::filesystem::remove_all;
using boost::filesystem::exists;

static void create_alignment(const string &name, const int position,
                             const string &cigar_string,
                             JewelerAlignment &al) {
    al.Position = position;
    al.Name = name;
    get_cigarop(cigar_string, al.CigarData);
    al.RefID = 0;
    al.AddTag("NM", "I", 1);
}

class SewingMachineTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        create_alignment("al1", 10, "10M", al1);
        create_alignment("al2", 20, "10M", al2);
        create_alignment("al3", 30, "10M", al3);
        create_alignment("al4", 40, "10M", al4);
        BamTools::RefData rf("chr1", 10);
        remove_all("test_data/test_sm/");
        this->initialize();
        add();
    }

    virtual void TearDown() {
        delete sm;
    }

    void initialize() {
        RefVector rv;
        rv.resize(10);
        string fd = "test_data/test_sm";
        sm = new SewingMachine;
        sm->load_zdb(fd);
        sm->references = rv;
    }

    void add() {
        sm->add_alignment(al1);
        sm->add_alignment(al2);
        sm->add_alignment(al2);
        sm->add_alignment(al3);
        al3.SetIsFirstMate(true);
        sm->add_alignment(al3);
        sm->add_alignment(al4);
        al3.SetIsFirstMate(true);
        sm->add_alignment(al4);
        sm->add_alignment(al4);
    }
    SewingMachine* sm;
    JewelerAlignment al1;
    JewelerAlignment al2;
    JewelerAlignment al3;
    JewelerAlignment al4;
};


TEST_F(SewingMachineTest, test_output_and_load) {
    ASSERT_EQ(false, sm->is_multiple_alignment("al1"));
    ASSERT_EQ(true, sm->is_multiple_alignment("al2"));
    ASSERT_EQ(false, sm->is_multiple_alignment("al3"));
    ASSERT_EQ(true, sm->is_multiple_alignment("al4"));
}


TEST_F(SewingMachineTest, test_output_save_and_load) {
    delete sm;
    this->initialize();
    ASSERT_EQ(false, sm->is_multiple_alignment("al1"));
    ASSERT_EQ(true, sm->is_multiple_alignment("al2"));
    ASSERT_EQ(false, sm->is_multiple_alignment("al3"));
    ASSERT_EQ(true, sm->is_multiple_alignment("al4"));
}
