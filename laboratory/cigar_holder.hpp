#include <iostream>
#include <string>
#include <api/BamReader.h>
#include <api/BamWriter.h>

using namespace BamTools;
using namespace std;
class JewelerAlignment;

int output_cigar_data(FILE * output_file, JewelerAlignment &al);
string get_cigar_string(const JewelerAlignment &al);
string get_cigar_string(const std::vector< CigarOp > &cigar_data);
void get_cigarop(const string &cigar_string, vector<CigarOp> &cigar_data);
void cigar_trim(JewelerAlignment &al);
