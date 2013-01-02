#ifndef _FASTA_H_
#define _FASTA_H_

#include <string>
#include <vector> 
using namespace std;

class seq_read {
public:

	string name;
	string seq;
	size_t size;
	seq_read(string name, string seq);
};

bool operator < (const seq_read &a, const seq_read &b);

/* 
   all existing data in the reads vector will be removed.
 */
int load_fasta_file(string fasta_filename, vector<seq_read *> &reads);

void write_fasta_file(const string filename,
                      const string &name,
                      const string &seq);
void write_fasta_file(FILE *,
                      const string &name,
                      const string &seq);

#endif /* _FASTA_H_ */
