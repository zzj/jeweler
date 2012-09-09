#ifndef _ZLEVELDB_
#define _ZLEVELDB_

#include "leveldb/db.h"
#include "leveldb/cache.h"
#include "common.hpp"
#include "proto/jeweler.pb.h"
#include <string>
#include <memory>
using namespace std;

class ZLevelDB {
public:
    ZLevelDB(string &db_folder) {
        this->db_folder = db_folder;
        initializ(db_folder);
    }

    void clear() {
        delete this->db;
        system((string("rm -rf ") + db_folder).c_str());
        system((string("mkdir -p ") + db_folder).c_str());
        // initialize again
        initializ(this->db_folder);
    }

    ~ZLevelDB() {
        delete this->db;
    }

    template<typename T>
    shared_ptr<T> get( const string &key);

    template<typename T>
    int set(const string &key, T *data);

    void del(string &key) {
        this->db->Delete(leveldb::WriteOptions(), key);
    }
private:
    leveldb::DB *db;
    string db_folder;

    void initializ(string &db_folder) {
        system((string("mkdir -p ") + db_folder).c_str());
        leveldb::Options options;
        options.create_if_missing = true;
        options.block_cache = leveldb::NewLRUCache(1000 * 1048576);
        leveldb::Status status = leveldb::DB::Open(options, db_folder, &this->db);
        assert(status.ok());
    }
};

template<typename T>
shared_ptr<T> ZLevelDB::get(const string &key) {
    string value;
    leveldb::Status s = this->db->Get(leveldb::ReadOptions(), key, &value);
    if (!s.ok()) return shared_ptr<T>();
    shared_ptr<T> ret(new T);
    if (!ret->ParseFromString(value)) {
        return shared_ptr<T>();
    }
    return ret;
}

template<typename T>
int ZLevelDB::set(const string &key, T *data) {
    string value;
    data->SerializeToString(&value);
    leveldb::Status s = this->db->Put(leveldb::WriteOptions(), key, value);
    if (!s.ok()) return -1;
    else return 0;
}

class ZMegaFile {
public:
    ZMegaFile(string &db_folder) : file_meta_db(db_folder) {
        this->file_name = db_folder + "/data.megafile";
        stream.open(file_name, ios::out | ios::binary | ios::in);
    }

    void clear() {
        stream.close();
        this->file_meta_db.clear();
        // make sure this function is called after db clear,
        // because db clear will clear all files under the folder.
        stream.open(file_name, ios::out | ios::binary | ios::in | ios::trunc);
    }

    template<typename T>
    shared_ptr<T> get(const string &key){
        shared_ptr<Jeweler::ZMegaFilePosition> zfp = \
            this->file_meta_db.get<Jeweler::ZMegaFilePosition>(key);
        if (zfp.get() == NULL) return shared_ptr<T>();
        shared_ptr<T> ret(new T);
        stream.seekg(zfp->position());
        load_protobuf_data(&stream, ret.get());
        return ret;
    }

    template<typename T>
    void append(const string &key, T*data){
        shared_ptr<Jeweler::ZMegaFilePosition> zfp(new Jeweler::ZMegaFilePosition());
        stream.seekg(0, ios_base::end);
        zfp->set_position(stream.tellg());
        this->file_meta_db.set(key, zfp.get());
        write_protobuf_data(&stream, data);
    }

    ~ZMegaFile() {
        stream.close();
    }

private:
    ZLevelDB file_meta_db;
    fstream stream;
    string file_name;
};

#endif
