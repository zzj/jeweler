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
        fprintf(stderr, "Failed to parse address book.\n");
        exit(1);
    }
    input.close();
}

#endif /* _COMMON_H_ */
