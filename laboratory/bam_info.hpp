#ifndef BAM_INFO
#define BAM_INFO

#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "main/common.hpp"
#include <map>
using namespace BamTools;
using namespace std;

class BamInfo{
public:
	// references
	SamHeader sam_header;
	RefVector references;
	virtual void initialize(BamReader &reader);

};

#endif 
