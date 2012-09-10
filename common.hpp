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
#include "jeweler_alignment.hpp"

using std::string;
using std::map;
using std::fstream;
using std::ios;

#define MAXLINE 10000

using namespace BamTools;

FILE * file_open(const char *name, const char * mode);

char * trim(char *);

void open_bam(BamReader &bam_reader, string bam_file);

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
    char * buffer = new char[m + 1];
    file->read(buffer, m);
    buffer[m] = '\0';
    // must use assign because '\0' might be part of the message.
    // TODO: add unit test.
    string raw; 
    raw.assign(buffer, m);
    if (!data->ParseFromString(raw)) {
        fprintf(stderr, "Failed to parse the buffer. [fstream version]");
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
    T *newdata = new T();
    return 0;
}

template<class T>
int map_add_count(map<T, int> &m, const T& key) {
    if (m.find(key) != m.end()) {
        m[key] ++;
    }
    else {
        m[key] = 0;
    }
    return m[key];
}

#endif /* _COMMON_H_ */
