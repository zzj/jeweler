#ifndef __GTF_INFO_H__
#define __GTF_INFO_H__
#include <vector>
#include <string>
using namespace std;

class gtf_info {
public:
    string type;
    unsigned int start, end;
    string chr;
    string transcript_id;
    string gene_id;
    gtf_info(vector<string> &strs);
};

#endif
