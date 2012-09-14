#ifndef _GTF_H_
#define _GTF_H_
#include <vector>
#include <string>
#include <Fasta.h>

using namespace std;

class Transcript;
void load_gtf_file(const string &gtf_filename,
                   vector<Transcript *> &maternal_transcripts,
                   vector<Transcript *> &paternal_transcripts,
                   FastaReference * maternal_ref,
                   FastaReference * paternal_ref);


#endif /* _GTF_H_ */
