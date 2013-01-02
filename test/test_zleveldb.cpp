#include "gtest/gtest.h"
#include "proto/jeweler.pb.h"
#include "main/zleveldb.hpp"

class ZLevelDBTest: public ::testing::Test {
protected:
    virtual void SetUp() {
        this->test_path = "test_data/test_zleveldb";
        this->zldb = new ZLevelDB(this->test_path);
    }

    virtual void TearDown() {
        this->zldb->clear();
        delete this->zldb;
    }

    void initialize_data() {
        Jeweler::ZMegaFilePosition *zfp = new Jeweler::ZMegaFilePosition;
        zfp->set_position(123);
        zldb->set("ABC", zfp);
        zfp->set_position(1234);
        zldb->set("ABCD", zfp);
    }
    
    ZLevelDB *zldb;
    string test_path;
};


TEST_F(ZLevelDBTest, test_set_and_get) {
    this->initialize_data();
    ASSERT_EQ(123, zldb->get<Jeweler::ZMegaFilePosition>("ABC")->position());
    ASSERT_EQ(1234, zldb->get<Jeweler::ZMegaFilePosition>("ABCD")->position());
}

TEST_F(ZLevelDBTest, test_set_save_and_get) {
    this->initialize_data();
    Jeweler::ZMegaFilePosition *zfp = new Jeweler::ZMegaFilePosition;
    zfp->set_position(1);
    ASSERT_EQ(0, zldb->set("ABC", zfp));
    ASSERT_EQ(1, zldb->get<Jeweler::ZMegaFilePosition>("ABC")->position());
    delete this->zldb;
    this->SetUp();
    ASSERT_EQ(1, zldb->get<Jeweler::ZMegaFilePosition>("ABC")->position());
    ASSERT_EQ(1234, zldb->get<Jeweler::ZMegaFilePosition>("ABCD")->position());
    ASSERT_EQ(NULL, zldb->get<Jeweler::ZMegaFilePosition>("ABCDE").get());
}

class ZMegaFileTest: public ::testing::Test {
protected:
    virtual void SetUp() {
        this->test_path = "test_data/test_zmegafile";
        this->zmf = new ZMegaFile(this->test_path);
    }

    virtual void TearDown() {
        this->zmf->clear();
        delete this->zmf;
    }

    void initialize_data() {
        Jeweler::ZMegaFilePosition *zfp = new Jeweler::ZMegaFilePosition;
        zfp->set_position(123);
        this->zmf->append("ABC", zfp);
        zfp->set_position(1234);
        this->zmf->append("ABCD", zfp);
    }

    void check() {
        ASSERT_TRUE(NULL != zmf->get<Jeweler::ZMegaFilePosition>("ABC").get());
        ASSERT_EQ(123, zmf->get<Jeweler::ZMegaFilePosition>("ABC")->position());
        ASSERT_TRUE(NULL != zmf->get<Jeweler::ZMegaFilePosition>("ABCD").get());
        ASSERT_EQ(1234, zmf->get<Jeweler::ZMegaFilePosition>("ABCD")->position());
        ASSERT_EQ(NULL, zmf->get<Jeweler::ZMegaFilePosition>("ABCDE").get());
    }

    ZMegaFile *zmf;
    string test_path;
};

TEST_F(ZMegaFileTest, test_append_and_get) {
    this->initialize_data();
    this->check();
}

TEST_F(ZMegaFileTest, test_append_save_and_get) {
    this->initialize_data();
    delete this->zmf;
    this->SetUp();
    this->check();
}

TEST_F(ZMegaFileTest, test_clear) {
    this->zmf->clear();
    this->zmf->stream.seekg(0, ios_base::end);
    ASSERT_EQ(0, this->zmf->stream.tellg());
}
