#ifndef _COMMON_H_
#define _COMMON_H_
#include <cstdio>

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <fstream>
#include <google/protobuf/io/coded_stream.h>
#include <memory>
using namespace google::protobuf::io;
using namespace std;
#define MAXLINE 10000

using namespace BamTools;

FILE * file_open(const char *name, const char * mode);

char * trim(char *);

void open_bam(BamReader &bam_reader, string bam_file);

double safe_rate(int a, int b);

#define _PROXY_(CLS) \
template <class T>                                              \
class read_only {                                                   \
    friend class CLS;                                           \
private:                                                        \
    T data;                                                     \
    T operator=(const T& arg) { data = arg; return data; }      \
    T operator+=(const T& arg) { data += arg; return data; }      \
    operator const T& () const { return data; }                \
    bool operator==(const T& arg) { return data == arg; }       \
    bool operator<=(const T& arg) { return data <= arg; }       \
    bool operator>=(const T& arg) { return data >= arg; }       \
    bool operator<(const T& arg) { return data < arg; }       \
    bool operator>(const T& arg) { return data > arg; }       \
    bool operator!=(const T& arg) { return data != arg; }       \
public:                                                         \
    const T& operator () () const { return data; }              \
};                                                              \

#define _PROXY_CLS_(CLS) \
template <class T>                                              \
class read_only_cls: private T {                                             \
    friend class CLS;                                           \
public:                                                         \
    const T& operator () () const { return *this; }              \
};                                                              \

template<class T>
void dump_protobuf_data(string file_name, T *data) {
    fstream out(file_name,
                ios::out | ios::binary | ios::trunc);
    data->SerializeToOstream(&out);
    out.close();
}

template<class T>
void load_protobuf_data(string file_name, T *data) {
    fstream input(file_name,
                ios::in | ios::binary);
    if (!data->ParseFromIstream(&input)) {
        fprintf(stderr, "Failed to parse the buffer. [string version]. FILE %s\n",
                file_name.c_str());
        exit(1);
    }
    input.close();
}

template<class T>
int load_protobuf_data(fstream *file, T *data) {
    int m;
    file->read(reinterpret_cast < char * > (&m), sizeof(m));
    if (file->eof()) return -1;
    ::google::protobuf::uint8 * buffer = new ::google::protobuf::uint8[m + 1];
    file->read((char *) buffer, m);
    // must use assign because '\0' might be part of the message.
    // TODO: add unit test.
    unique_ptr<CodedInputStream> input(new CodedInputStream(buffer, m));
    input->SetTotalBytesLimit(1024 * 1024 * 1024, 1024 * 1024 * 1024);
    if (!data->ParseFromCodedStream(input.get())) {
        fprintf(stderr, "Failed to parse the buffer. [fstream version]");
        delete buffer;
        return -1;
    }
    delete buffer;
    return 0;
}

template<class T>
int write_protobuf_data(fstream *file, T *data) {
    string seq;
    data->SerializeToString(&seq);
    int m = seq.size();
    file->write(reinterpret_cast<char *> (&m), sizeof(m));
    file->write(seq.c_str(), seq.size());
    return 0;
}

template<class T>
int map_add_count(map<T, int> &m, const T& key) {
    if (m.find(key) != m.end()) {
        m[key] ++;
    }
    else {
        m[key] = 1;
    }
    return m[key];
}

template<class T, class U>
bool map_add_default(map<T, U> &m, const T& key, const U& value) {
    if (m.find(key) != m.end()) {
        return false;
    }
    else {
        m[key] = value;
        return true;
    }
}

template<class T, class U>
U map_get_default(const map<T, U> &m, const T& key, const U& value) {
    typename map<T,U>::const_iterator it;
    if ((it = m.find(key)) != m.end()) {
        return it->second;
    }
    else {
        return value;
    }
}

template<class T>
vector<T *> duplicate_vector(vector<T *> in) {
    vector<T *> ret(in.size());
    size_t i;
    for (i = 0; i < in.size(); i ++) {
        ret[i] = new T(*in[i]);
    }
    return ret;
}

template<class T, class U>
U map_value_sum(class map<T, U> &m, U init) {
    for (auto i = m.begin(); i != m.end(); i ++) {
        init += i->second;
    }
    return init;
}


#endif /* _COMMON_H_ */
