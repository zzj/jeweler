#include "cigar_holder.hpp"
#include <sstream>
using namespace std;
int output_cigar_data(FILE * output_file, JewelerAlignment &al) {
	std::vector< CigarOp > &cigar_data = al.CigarData;
	for (auto i=cigar_data.begin(); i!=cigar_data.end(); i++) {
		fprintf(output_file,"%d%c",i->Length,i->Type);
	}
	fprintf(output_file,"\n");
	return 0;
}

string  get_cigar_string(const JewelerAlignment &al) {
	const std::vector< CigarOp > &cigar_data = al.CigarData;
	ostringstream oss;
	for (auto i=cigar_data.begin(); i!=cigar_data.end(); i++) {
		oss << i->Length << i->Type;
	}	
	return oss.str();
}
