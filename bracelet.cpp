#include "bracelet.hpp"
#include "common.hpp"
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

	reads.resize(jeweler_info->gene_id2transcripts.size());
	reads_index.resize(jeweler_info->gene_id2transcripts.size());
	results.resize(jeweler_info->gene_id2transcripts.size());
	related_transcripts.resize(jeweler_info->gene_id2transcripts.size());
    for (auto i = jeweler_info->gene_id2transcripts.begin();
         i != jeweler_info->gene_id2transcripts.end();
         i++) {
        string gene_id = i->first;
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
        map_add_count(coverage, pos);
    }
    return ;
}


double get_coverage_rate(const Jeweler::EarringsData::Read &origin,
                         const map<int, int> &coverage) {
    if (origin.genome_position_size() == 0) return 0;
    int cov = 0;
    for (int i = 0 ; i < origin.genome_position_size(); i ++) {
        const int pos = origin.genome_position(i);
        if (map_get_default(coverage, pos, 0) > 0) {
            cov ++;
        }
    }
    return double(cov) / origin.genome_position_size();
}

void add_coverage_details(const Jeweler::EarringsData::Read &origin,
                          const Jeweler::EarringsData::Read &target,
                          map<int, map<int, int> > &genome_position_map) {
    // genome_position_map: first is the original genome position,
    // second is another map, which is the target_genome_position and times.
    map<int, int> original_read2genome;
    assert(origin.genome_position_size() == origin.read_position_size());

    for (int i = 0; i < origin.genome_position_size(); i ++){
        original_read2genome[origin.read_position(i)] = origin.genome_position(i);
    }
    assert(target.genome_position_size() == target.read_position_size());
    for (int i = 0; i < target.genome_position_size(); i ++){
        if (original_read2genome.find(target.read_position(i)) != 
            original_read2genome.end()) {
            map_add_default(genome_position_map,
                            target.read_position(i),
                            map<int, int>());
            map_add_count(genome_position_map[target.read_position(i)],
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

    for (size_t i = 0; i < ed->read_size(); i++) {
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

    for (size_t i = 0; i < ed->read_size(); i++) {
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
