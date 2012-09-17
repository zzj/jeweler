#include "bracelet.hpp"
#include "common.hpp"
#include "constants.hpp"
#include "jeweler_info.hpp"
#include <boost/filesystem.hpp>
#include "zleveldb.hpp"

using boost::filesystem::path;
using boost::filesystem::create_directory;

Bracelet::Bracelet(JewelerInfo * jeweler_info) :
    zmf(new ZMegaFile(jeweler_info->result_folder)) {
	int id = 0;
    this->num_single_reads = 0;
    this->num_fixed_reads = 0;

	this->jeweler_info = jeweler_info;
	fprintf(stdout, "bracelet initializing ...\n");

	reads.resize(jeweler_info->get_num_genes());
	reads_index.resize(jeweler_info->get_num_genes());
	results.resize(jeweler_info->get_num_genes());
	related_transcripts.resize(jeweler_info->get_num_genes());
    for (auto i = jeweler_info->gene_id.begin();
         i != jeweler_info->gene_id.end();
         i++) {
        string gene_id = (*i);
        shared_ptr<Jeweler::EarringsData> ed = 
            this->zmf->get<Jeweler::EarringsData>(gene_id);
        if (ed.get() == NULL) continue;
        for (int j = 0; j < ed->read_size(); j++) {
            if (ed->read(j).is_multiple_alignment()) {
                reads[id].push_back(ed->read(j).name());
                reads_index[id][ed->read(j).name()] = j;
            }
        }
        sort(reads[id].begin(),reads[id].end());
        id ++;
    }
    this->data = new Jeweler::BraceletData();
}

int Bracelet::intersect(vector<string> &a, vector<string>&b) {
    size_t i , j , r;
    i = 0, j = 0;
    r = 0;
    while(i != a.size() && j != b.size()) {
        if (a[i]<b[j]) {i++; }
        else if (a[i]>b[j]) {j++; }
        else {i++,j++,r++;}
    }
    return r;
}


int Bracelet::analyze() {
    fprintf(stdout, "bracelet analyzing ...\n");
    for (size_t i=0;i<reads.size(); i++) {
        for (size_t j=i+1;j<reads.size();j++) {
            int r = intersect(reads[i], reads[j]);
            if (r > 0) {
                results[i].push_back(r);
                related_transcripts[i].push_back(j);
                results[j].push_back(r);
                related_transcripts[j].push_back(i);
            }
        }
    }
    return 0;
}

void add_coverage(const Jeweler::EarringsData::Read &origin,
                  map<int, int> &coverage) {
    for (int i = 0 ; i < origin.genome_position_size(); i ++) {
        int pos = origin.genome_position(i);
        if (pos != NOT_FOUND)
            map_add_count(coverage, pos);
    }
    return ;
}


double get_coverage_rate(const Jeweler::EarringsData::Read &origin,
                         const map<int, int> &coverage) {
    if (origin.genome_position_size() == 0) return 0;
    int cov = 0;
    int total = 0;
    for (int i = 0 ; i < origin.genome_position_size(); i ++) {
        const int pos = origin.genome_position(i);
        if (pos == NOT_FOUND) continue;
        if (map_get_default(coverage, pos, 0) > 0) {
            cov ++;
        }
        total ++;
    }
    return double(cov) / total;
}


// This is arbitrarily implemented.

//                ~~~~~~  Ambiguous region.
// Origin1: ------------
// Origin2:       ************
// HEAD:          ^
// TAIL:                ^
//          ~~~~~~                  ~~~~~
// Target1: ************
// Target2:                  ____________

int get_origin_read_position(const Jeweler::EarringsData::Read &origin,
                             const Jeweler::EarringsData::Read &target,
                             int i) {
    int origin_overlap_head, origin_overlap_tail, target_overlap_head, target_overlap_tail;
    origin_overlap_tail = origin.glue_position();
    origin_overlap_head = origin.genome_position_size() - origin.tail_length();
    target_overlap_tail = target.glue_position();
    target_overlap_head = target.genome_position_size() - target.tail_length();
    if (origin.is_second_truncated() == target.is_second_truncated()) {
        if (i < target.head_length()) {
            return i;
        }
        else {
            return (i - target_overlap_head) + origin_overlap_head;
        }
    }
    else {
        if (i < target.head_length()) {
            return i + origin_overlap_head;
        }
        else {
            return i - target.head_length();
        }
     }
}
                          

void add_coverage_details(const Jeweler::EarringsData::Read &origin,
                          const Jeweler::EarringsData::Read &target,
                          map<int, map<int, int> > &genome_position_map) {
    // genome_position_map: first is the original genome position,
    // second is another map, which is the target_genome_position and times.
    map<int, int> original_read2genome;

    for (int i = 0; i < origin.genome_position_size(); i ++){
        if (origin.genome_position(i) != NOT_FOUND)
            original_read2genome[i] = origin.genome_position(i);
    }
    for (int i = 0; i < target.genome_position_size(); i ++){
        if (target.genome_position(i) == NOT_FOUND)
            continue;
        int origin_i = get_origin_read_position(origin, target, i);
        if (original_read2genome.find(origin_i) != 
            original_read2genome.end()) {
            map_add_default(genome_position_map,
                            i,
                            map<int, int>());
            map_add_count(genome_position_map[i],
                          target.genome_position(i));
        }
    }
    return ;
}


void Bracelet::dump_shared_pileup(Jeweler::BraceletData::RelatedTranscript * rt,
                                  int original_id,
                                  int target_id) {

    map<int, map<int, int> > coverage_details;
    map<int, int> coverage;
    shared_ptr<Jeweler::EarringsData> ed = 
        this->zmf->get<Jeweler::EarringsData>(jeweler_info->gene_id[original_id]);

    if (ed.get() == NULL)
        return ;

    shared_ptr<Jeweler::EarringsData> target_ed = 
        this->zmf->get<Jeweler::EarringsData>(jeweler_info->gene_id[target_id]);

    if (target_ed.get() == NULL)
        return ;

    for (int i = 0; i < ed->read_size(); i++) {
        string name = ed->read(i).name();
        if (reads_index[original_id].find(name)==reads_index[original_id].end() ||
            reads_index[target_id].find(name)==reads_index[target_id].end()) {
            continue;
        }
        add_coverage(ed->read(i), coverage);
        add_coverage_details(ed->read(i),
                             target_ed->read(reads_index[target_id][name]),
                             coverage_details);
    }

    for (int i = 0; i < ed->read_size(); i++) {
        if (ed->read(i).is_multiple_alignment()) {
            continue;
        }
        else {
            double rate = get_coverage_rate(ed->read(i), coverage);
            num_single_reads ++;
            if (rate > 0.9) {
                num_fixed_reads ++;
            }
        }
    }

    for (auto i = coverage.begin(); i != coverage.end(); i ++) {
        Jeweler::BraceletData::RelatedTranscript::Coverage *c = rt->add_coverage();
        c->set_position(i->first);
        c->set_total_coverage(i->second);
        for (auto j = coverage_details[i->first].begin();
             j != coverage_details[i->first].end();
             j ++) {
            c->add_target_position(j->first);
            c->add_shared_coverage(j->second);
        }
    }
}


int Bracelet::dump(fstream *fd, string root) {
    fprintf(stdout, "bracelet dumping %s...\n", root.c_str());
    for (size_t i = 0; i < reads.size() ;i++) {
        if (results[i].size()==0) continue;
        unique_ptr<Jeweler::BraceletData> t(new Jeweler::BraceletData());
        t->set_name(jeweler_info->gene_id[i]);
        for (size_t j = 0; j < results[i].size(); j ++) {
            Jeweler::BraceletData::RelatedTranscript *rt = t->add_related_transcript();
            rt->set_name(jeweler_info->gene_id[related_transcripts[i][j]]);
            rt->set_num_shared_read(results[i][j]);
            dump_shared_pileup(rt, i, related_transcripts[i][j]);
		}
        write_protobuf_data(fd, t.get());
	}
    fprintf(stdout, "%d/%d", this->num_fixed_reads, this->num_single_reads);
	return 0;
}
