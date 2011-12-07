#ifndef _GTF_H_
#define _GTF_H_
#include "transcript.hpp"
#include <boost/algorithm/string.hpp>
int load_gtf_file(string gtf_filename, vector<Transcript *> &transcripts);

#endif /* _GTF_H_ */
