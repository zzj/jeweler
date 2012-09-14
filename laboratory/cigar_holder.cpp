#include "cigar_holder.hpp"
#include <sstream>
#include "../common.hpp"
#include "../jeweler_alignment.hpp"
using namespace std;
int output_cigar_data(FILE * output_file, JewelerAlignment &al) {
	std::vector< CigarOp > &cigar_data = al.CigarData;
	for (auto i=cigar_data.begin(); i!=cigar_data.end(); i++) {
		fprintf(output_file,"%d%c",i->Length,i->Type);
	}
	fprintf(output_file,"\n");
	return 0;
}

string  get_cigar_string(const std::vector< CigarOp > &cigar_data) {
	ostringstream oss;
	for (auto i=cigar_data.begin(); i!=cigar_data.end(); i++) {
		oss << i->Length << i->Type;
	}
	return oss.str();
}

string  get_cigar_string(const JewelerAlignment &al) {
	return get_cigar_string(al.CigarData);
}

void get_cigarop(const string &cigar_string, vector<CigarOp> &cigar_data) {
    istringstream iss;
    char Type;
    int Length;
    iss.str(cigar_string);
    cigar_data.clear();
    while(iss >> Length >> Type) {
        CigarOp co;
        co.Type = Type;
        co.Length = Length;
        cigar_data.push_back(co);
    }
    return;
}

void cigar_trim(JewelerAlignment &al) {

	vector<CigarOp>& cigar_data = al.CigarData;
	vector<CigarOp> new_cigar_data;

	auto i = cigar_data.begin();
	// trim the beginning 'N'
	while( i != cigar_data.end() && i->Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
        al.Position += i->Length;
        i ++;
    }

	while ( i != cigar_data.end() ) {
		auto j = i;
		j++;
		if (j != cigar_data.end() ) {
			if (i->Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
				while( j !=cigar_data.end() &&
					   j->Type == Constants::BAM_CIGAR_REFSKIP_CHAR
					   ) {
					i->Length += j->Length;
					j++;
				}
			}
		}
		if ((i->Type != Constants::BAM_CIGAR_REFSKIP_CHAR || j != cigar_data.end() ))
			new_cigar_data.push_back(*i);
		i = j;
	}
	al.CigarData = new_cigar_data;
	return ;
}
