#ifndef _GTF_H_
#define _GTF_H_
#include "transcript.hpp"
#include <boost/algorithm/string.hpp>
void load_gtf_file(const string &gtf_filename, vector<Transcript *> &transcripts);

#endif /* _GTF_H_ */
