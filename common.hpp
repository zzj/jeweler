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
#include "jeweler_alignment.hpp"


#define MAXLINE 10000

using namespace BamTools;

FILE * file_open( const char *name, const char * mode);

char * trim(char *);

void open_bam(BamReader &bam_reader, string bam_file);

#endif /* _COMMON_H_ */
