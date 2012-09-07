#ifndef _ZLEVELDB_
#define _ZLEVELDB_

#include "leveldb/db.h"
#include "leveldb/cache.h"
#include <string>
#include <memory>
using namespace std;
class ZLevelDB {
public:
    ZLevelDB(string &db_folder) {
        initializ(db_folder);
    }

    void initializ(string &db_folder) {
        system((string("mkdir -p ") + db_folder).c_str());
        leveldb::Options options;
        options.create_if_missing = true;
        options.block_cache = leveldb::NewLRUCache(1000 * 3048576);
        leveldb::Status status = leveldb::DB::Open(options, db_folder, &this->db);
        assert(status.ok());
    }

    void clear(string &db_folder) {
        delete this->db;
        system((string("rm -rf ") + db_folder).c_str());
        initializ(db_folder);
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

    leveldb::DB *db;
};

template<typename T>
shared_ptr<T> ZLevelDB::get(const string &key) {
    string value;
    leveldb::Status s = this->db->Get(leveldb::ReadOptions(), key, &value);
    if (!s.ok()) return shared_ptr<T>(new T);
    shared_ptr<T> ret(new T);
    if (!ret->ParseFromString(value)) {
        return shared_ptr<T>(new T);
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
#endif
