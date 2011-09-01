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


#define MAXLINE 10000

FILE * file_open( const char *name, const char * mode);

char * trim(char *);
#endif /* _COMMON_H_ */
