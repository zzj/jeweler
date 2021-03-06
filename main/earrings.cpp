#include "earrings.hpp"
#include "jeweler_info.hpp"
#include "proto/jeweler.pb.h"
#include "laboratory/sewing_machine.hpp"
#include "aligner.hpp"
#include "transcript_mismatcher.hpp"
#include "pileup.plot.hpp"
#include "alignment_glue.hpp"
#include "graph/graph.hpp"
#include "zleveldb.hpp"
#include "constants.hpp"
#include <fstream>
using namespace std;

Earrings::Earrings(JewelerInfo *jeweler_info,
                   string gene_id,
                   SewingMachine *sm,
                   ZMegaFile *file) {
    size_t i;

    this->sm = sm;

    this->jeweler_info = jeweler_info;
    this->gene_id = gene_id;
    this->result_folder = jeweler_info->result_folder + "/" + gene_id + "/";

    if (!boost::filesystem::create_directory(this->result_folder))
        fprintf(stderr, "cannot create folder %s", this->result_folder.c_str());
    this->zmf = file;

    this->maternal_transcripts =
        this->jeweler_info->get_maternal_transcripts(this->gene_id);
    this->paternal_transcripts =
        this->jeweler_info->get_paternal_transcripts(this->gene_id);

    this->mismatcher = new TranscriptMismatcher();
    for (i = 0; i < paternal_transcripts.size(); i++) {
        mismatcher->add_transcript(paternal_transcripts[i], TRANSCRIPT_PATERNAL);
        mismatcher->add_transcript(maternal_transcripts[i], TRANSCRIPT_MATERNAL);
    }
    mismatcher->initialize();

    this->chr = paternal_transcripts[0]->chr();
    this->ref_id = jeweler_info->get_refID(chr);
    if (ref_id == NOT_FOUND) {
        fprintf(stdout, "Cannot find reference id for %s\n", this->chr.c_str());
    }

    for (i=0;i<paternal_transcripts.size();i++) {
        Transcript *p=paternal_transcripts[i];

        if (i == 0) {
            left_pos = p->start();
            right_pos = p->end();
        }
        else {
            left_pos = min(left_pos, p->start());
            right_pos = max(right_pos, p->end());
        }
    }

    // load bamalignment by bamtools.
    load_read_data();

    // align reads back to the transcripts
    align_reads();
    // build graph;
    //build_graph();
    //test_allele_specific_transcript();
    // have sewingmachine class, output maltiple aligment info
    dump_data();
}

Earrings::~Earrings() {
    size_t i=0;
    test_memory_leak();
    for (i = 0; i < this->maternal_transcripts.size(); i ++) {
        // hack, because transcripts objects are becoming larger
        // during the running time.
        // need to remove some resources.
        delete this->maternal_transcripts[i];
        delete this->paternal_transcripts[i];
    }
    for (i = 0; i < bam_reads.size(); i++) {
        delete bam_reads[i];
    }
    for (i = 0; i < unaligned_reads.size(); i++) {
        delete unaligned_reads[i];
    }
    for (i = 0; i < noused_reads.size(); i++) {
        delete noused_reads[i];
    }
    delete mismatcher;
}



void Earrings::test_memory_leak() {
    if (num_total_reads !=
        bam_reads.size() + unaligned_reads.size() + noused_reads.size())
        fprintf (stderr, "WARNING: Memory Leaking: total reads %zu\t record reads%zu\n",
                 num_total_reads,
                 bam_reads.size() + unaligned_reads.size() + noused_reads.size());
}

int Earrings::load_read_data() {
    num_total_reads = 0;

    if (! jeweler_info->bam_reader.SetRegion(ref_id, left_pos, ref_id, right_pos)) {
        fprintf(stdout, "%s\n", jeweler_info->bam_reader.GetErrorString().c_str());
    }

    JewelerAlignment *al=new JewelerAlignment();

    while(jeweler_info->bam_reader.GetNextAlignment(*al)) {
        al->jeweler_initialize(this->sm);
        bam_reads.push_back(al);
        al=new JewelerAlignment();
    }
    delete al; // delete the last unused one
    num_total_reads = bam_reads.size();

    fprintf(stdout, "totally %zu bamalignments are loaded\n", num_total_reads);

    return 0;
}

int Earrings::study_compatible_reads() {
    Transcript::tolerate = 0;
    for (int t = 0; t < 10; t++) {
        Transcript::tolerate=t;
        //get_compatible_reads();
    }
    return 0;
}


void Earrings::get_compatible_reads(vector<set<JewelerAlignment *> > &read_lists) {
    size_t i,j;
    size_t num_unaligned_reads = 0;
    bool is_compatible=false;
    bool is_fixed;
    read_lists.resize(maternal_transcripts.size());
    compatible_reads.clear();
    for (i = 0; i < bam_reads.size(); i++) {
        int min_penalty =  100000;
        vector<int> penalties;
        penalties.resize(maternal_transcripts.size(),100000);
        is_compatible=false;

        for (j = 0; j < maternal_transcripts.size(); j++) {
            // the maternal_transcript should be the same with the
            // paternal_transcript

            if (maternal_transcripts[j]->is_compatible(bam_reads[i],
                                                       Transcript::tolerate)) {
                maternal_transcripts[j]->get_overlapped_alignment(bam_reads[i],
                                                                  penalties[j]);
                is_compatible=true;
                if (penalties[j] < min_penalty) {
                    min_penalty = penalties[j];
                }
            }
        }
        is_fixed=false;

        if (!is_compatible  || min_penalty > Transcript::tolerate) {
            unaligned_reads.push_back(bam_reads[i]);
            num_unaligned_reads ++;
        }
        else {
            compatible_reads.push_back(bam_reads[i]);
            for (j = 0; j < maternal_transcripts.size(); j++) {
                if (penalties[j] == min_penalty) {
                    if (min_penalty != 0 && !is_fixed) {
                        maternal_transcripts[j]->get_overlapped_alignment(bam_reads[i],
                                                                          penalties[j],
                                                                          true);
                        is_fixed = true;
                    }
                    else {
                        maternal_transcripts[j]->get_overlapped_alignment(bam_reads[i],
                                                                          penalties[j],
                                                                          false);
                        // oh shit, the two transcripts are have
                        // different ways to make them compatible
                        if (penalties[j] > 0) {
                            continue;
                        }
                    }
                    read_lists[j].insert(bam_reads[i]);
                }
            }
        }
    }

    fprintf(stdout, "Tolerate %d\tTotal reads: %zu\t Compatible: %zu\t unaligned: %zu\n",
            Transcript::tolerate,
            bam_reads.size(),
            bam_reads.size()  - num_unaligned_reads,
            num_unaligned_reads);


}

void Earrings::align_reads() {
    size_t i,j,k;

    int total_alleles;
    bool is_compatible=false;
    set<JewelerAlignment *> cleared;
    vector<set<JewelerAlignment *> >noninfo;

    vector<JewelerAlignment *> new_bam_reads;
    vector<set<JewelerAlignment *> > read_lists;

    AlignmentGlue ag;
    Transcript::tolerate = 0;

    ag.glue(bam_reads, new_bam_reads, noused_reads);
    printf("%zu new paired reads, %zu unused reads\n",
           new_bam_reads.size(), noused_reads.size());
    bam_reads=new_bam_reads;

    get_compatible_reads(read_lists);

    bam_reads=compatible_reads;

    int num_compatible_reads = 0;
    noninfo.resize(maternal_transcripts.size());
    for (i = 0 ; i < bam_reads.size(); i++) {

        is_compatible = false;
        for (j = 0; j < maternal_transcripts.size(); j++) {
            if (read_lists[j].find(bam_reads[i]) == read_lists[j].end()) {
                continue;
            }
            // the maternal_transcript should be the same with the paternal_transcript
            if (maternal_transcripts[j]->is_compatible(bam_reads[i], Transcript::tolerate)) {
                is_compatible = true;

                // TODO: use a class to store the matched information
                ReadMatcher maternal_matcher;
                ReadMatcher paternal_matcher;

                maternal_transcripts[j]->match_alleles(bam_reads[i],
                                                       total_alleles,
                                                       &maternal_matcher);

                paternal_transcripts[j]->match_alleles(bam_reads[i],
                                                       total_alleles,
                                                       &paternal_matcher);

                //fprintf(stdout,"%d\n",total_alleles);
                int num_maternal_alleles =
                    maternal_matcher.allele_read_locations.size();
                int num_paternal_alleles =
                    paternal_matcher.allele_read_locations.size();

                if (total_alleles>0 &&
                    num_paternal_alleles != num_maternal_alleles) {

                    cleared.insert(bam_reads[i]);
                    if (num_maternal_alleles > num_paternal_alleles) {
                        maternal_transcripts[j]->register_allele_read(bam_reads[i],
                                                                      maternal_matcher);
                        mismatcher->add_mismatches(maternal_transcripts[j],
                                                   bam_reads[i],
                                                   &maternal_matcher);
                    }
                    else{
                        paternal_transcripts[j]->register_allele_read(bam_reads[i],
                                                                      paternal_matcher);
                        mismatcher->add_mismatches(paternal_transcripts[j],
                                                   bam_reads[i],
                                                   &paternal_matcher);
                    }
                }
                else {
                    maternal_transcripts[j]->register_read(bam_reads[i]);
                    paternal_transcripts[j]->register_read(bam_reads[i]);
                    mismatcher->add_mismatches(maternal_transcripts[j],
                                               bam_reads[i],
                                               &maternal_matcher);
                    noninfo[j].insert(bam_reads[i]);
                }
            }
        }

        if (is_compatible) {
            num_compatible_reads ++;
        }
        else {
            fprintf(stderr, "WARNING: Find an incompatible read.\n");
            for (k = 0; k < paternal_transcripts.size(); k ++) {
                paternal_transcripts[k]->output_segments();
                paternal_transcripts[k]->is_compatible(bam_reads[i], Transcript::tolerate,
                                                   /*debug output*/true);
            }
            output_bamalignment(bam_reads[i]);
        }
    }

    fprintf(stdout, "%d compatible reads in %zu reads\n",
            num_compatible_reads,
            bam_reads.size()
            );

    FILE *finfo = fopen(string(result_folder+"/"+gene_id+".landscape.plot.meta").c_str(),"w+");
    for (i=0;i<maternal_transcripts.size();i++) {
        FILE *foutput;
        auto filename = \
            result_folder + "/" + maternal_transcripts[i]->transcript_id() + \
            ".landscape.plot.info";
        foutput = fopen(filename.c_str(),"w+");
        if (foutput == NULL) {
            fprintf(stderr, "cannot open file at %s\n",
                    filename.c_str());
        }
        string name=maternal_transcripts[i]->transcript_id();
        PileupPlot lp(maternal_transcripts[i],
                      paternal_transcripts[i],
                      noninfo[i], multiple_reads);
        lp.generate_pileup_plot(finfo, foutput);
        fclose(foutput);
    }
    fclose(finfo);

    //fprintf(stdout,"Unaligned\tUncleared\tCleared\tNoninfo\tTotal\n");
    //fprintf(stdout,"%d\t%d\t%d\t%d\t%d\n",unaligned.size(),
    //uncleared.size(),
    //cleared.size(),noninfo.size(),bam_reads.size());

}

int Earrings::test_allele_specific_transcript() {
    size_t i,j;
    bool maternal_dominance=false, paternal_dominance=false;
    for (i = 0; i < maternal_transcripts.size(); i++) {

        vector<int>& maternal=maternal_transcripts[i]->num_info_reads_per_exon;
        vector<int>& paternal=paternal_transcripts[i]->num_info_reads_per_exon;
        for (j = 0; j < maternal.size(); j++) {
            //if (maternal[j]>(paternal[j]+1)*3) {
            if (maternal[j] > paternal[j]*10) {
                maternal_dominance=true;
            }
            //if (paternal[j]>(maternal[j]+1)*3) {
            if (paternal[j]>10* maternal[j]) {
                paternal_dominance=true;
            }
        }

    }
    return (maternal_dominance || paternal_dominance);

}


int Earrings::build_graph() {
    size_t i,j;
    Graph graph;
    fprintf(stdout,"%s\n",gene_id.c_str());
    FILE *foutput=fopen(string(result_folder+"/"+gene_id+".allele.specific.graph").c_str(),"w+");
    vector<Path> cufflink_records;
    for (i=0;i<maternal_transcripts.size();i++) {
        maternal_transcripts[i]->add_transcript_to_graph(&graph,cufflink_records);
        paternal_transcripts[i]->add_transcript_to_graph(&graph,cufflink_records);
    }
    graph.dump_graph(foutput);
    fclose(foutput);

    vector<Path> records;
    graph.get_all_paths(records);
    // get allele specific isoforms
    for (i = 0; i < records.size(); i++) {
        bool is_mirrored=false;
        for (j = 0; j < records.size(); j++) {
            if (records[i].is_mirrored(records[j])) {
                is_mirrored=true;
                break;
            }
        }
        if (!is_mirrored &&records[i].is_informative()) {
            //records[i].dump_path(stdout);
        }
    }

    for (i = 0; i < cufflink_records.size(); i++) {
        if (!cufflink_records[i].is_valid())
            cufflink_records[i].dump_path(stdout);
    }
    //get allele specific isoforms that is not from cufflinks

    FILE *foutput_path=fopen(string(result_folder+"/"+gene_id+".allele.specific.path").c_str(),"w+");
    for (j = 0; j < records.size();j++) {
        bool is_new=true;
        for (i = 0; i < cufflink_records.size(); i++) {
            if (cufflink_records[i].is_valid()) {
                if (cufflink_records[i].is_equal(records[j])) {
                    is_new=false;
                }
            }
        }
        if (is_new) {
            records[j].dump_path(foutput_path);
        }
    }
    fclose(foutput_path);
    FILE *foutput_info=fopen(string(result_folder+"/"+gene_id+".allele.specific.info").c_str(),"w+");
    if (test_allele_specific_transcript()) {
        fprintf(foutput_info,"%s","Yes");
    }
    else {
        fprintf(foutput_info,"%s","No");
    }
    fclose(foutput_info);
    return 0;
}

void Earrings::dump_data() {

    unique_ptr<Jeweler::EarringsData> ed(new Jeweler::EarringsData());
    ed->set_gene_id(this->gene_id);
    ed->set_chr(this->chr);
    ed->set_start_position(this->left_pos);
    ed->set_end_position(this->right_pos);
    ed->set_num_single_reads(this->single_reads.size());
    ed->set_num_multiple_reads(this->multiple_reads.size());
    this->mismatcher->dump(ed->mutable_mismatcher());
    this->dump_reads(ed.get());
    for (int i = 0; i < maternal_transcripts.size(); i++) {
        auto transcript_data = ed->add_transcript();
        transcript_data->set_transcript_id(maternal_transcripts[i]->transcript_id());
        transcript_data->set_maternal_seq(maternal_transcripts[i]->seq());
        transcript_data->set_paternal_seq(paternal_transcripts[i]->seq());
        for (int j = 0; j < maternal_transcripts[i]->exon_start.size(); j ++) {
            auto exon_data = transcript_data->add_exon();
            exon_data->set_start(maternal_transcripts[i]->exon_start[j]);
            exon_data->set_end(maternal_transcripts[i]->exon_end[j]);
        }
    }
    this->zmf->append(this->gene_id, ed.get());
    return ;
}

void Earrings::classify_reads() {

    single_reads.clear();
    multiple_reads.clear();

    for (size_t i = 0; i < bam_reads.size(); i ++) {
        if (sm->is_multiple_alignment(bam_reads[i]->Name)){
            multiple_reads.insert(bam_reads[i]);
        }
        else {
            single_reads.insert(bam_reads[i]);
        }
    }

    return ;
}

void Earrings::dump_reads(Jeweler::EarringsData *ed) {
    for (auto i = compatible_reads.begin(); i != compatible_reads.end(); i++) {
        (*i)->dump_data(ed->add_read());
    }
    return ;
}
