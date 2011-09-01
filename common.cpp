
#include "common.hpp"
using namespace std;

char *trim(char *str)
{
	char *end;

	// Trim leading space
	while(isspace(*str)) str++;

	if(*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace(*end)) end--;

	// Write new null terminator
	*(end+1) = 0;

	return str;
}

FILE * file_open( const char *name, const char * mode){
     FILE * ret=NULL;
     ret=fopen(name, mode);
     if (ret==NULL)  {
	  fprintf(stderr, "ERROR: can not open file %s  at %s:%d\n",
			     name, __FILE__, __LINE__);
	  exit(0);
     }
     return ret;
}
