#include <iostream>
#include <string> 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "../common.hpp"

using namespace BamTools;
using namespace std;

int output_cigar_data(FILE * output_file, BamAlignment &al);
string get_cigar_string(BamAlignment &al);
