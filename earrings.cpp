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
#include <fstream>
using namespace std;

Earrings::Earrings(JewelerInfo *jeweler_info,
				   string gene_id,
				   SewingMachine *sm,
                   ZMegaFile *file,
				   bool is_prepare = false) {
	size_t i;

	this->sm = sm;
	this->mismatcher = new TranscriptMismatcher();
	this->jeweler_info = jeweler_info;
	this->gene_id = gene_id;
	this->result_folder = jeweler_info->result_folder + "/" + gene_id + "/";
    this->zmf = file;

	// load maternal and paternal transcripts sequences.
	load_transcript_data(is_prepare);
	//if (is_prepare) return;

	this->chr = paternal_transcripts[0]->chr;
	this->ref_id = jeweler_info->get_refID(chr);
	if (ref_id == NOT_FOUND) {
		fprintf(stdout, "Cannot find reference id for %s\n", this->chr.c_str());
	}

	for (i=0;i<paternal_transcripts.size();i++) {
		Transcript *p=paternal_transcripts[i];

		if (i == 0) {
			left_pos = p->start;
			right_pos = p->end;
		}
		else {
			left_pos = min(left_pos, p->start);
			right_pos = max(right_pos, p->end);
		}
	}
	// load bamalignment by bamtools.
	load_read_data();
	if (sm!=NULL) {
		count_multiple_alignments(/* is after alignment*/ false);
	}

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
	for (i = 0; i < maternal_transcripts.size(); i++) {
		delete maternal_transcripts[i];
		delete paternal_transcripts[i];
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


// this function must be called after the align_reads
// otherwise the statistics may be incorrect.
void Earrings::count_multiple_alignments(bool is_after_aligned) {

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

    FILE *foutput_mamf = fopen(string(result_folder+"/"+gene_id+".mamf.meta").c_str(),"w+");
    fprintf(foutput_mamf, "%zu\t%zu\n", single_reads.size(), bam_reads.size());
    fprintf(stdout, "%zu\t%zu\n", single_reads.size(), bam_reads.size());
    fclose(foutput_mamf);
    if (is_after_aligned) {
        foutput_mamf=fopen(string(result_folder+"/"+gene_id+".mamf.multiple.reads").c_str(),"w+");
    }
    else {
        foutput_mamf=fopen(string(result_folder+"/"+gene_id+".mamf.before.multiple.reads").c_str(),"w+");
    }
    for (auto i = multiple_reads.begin(); i != multiple_reads.end(); i++) {
        fprintf(foutput_mamf, "%s\t%s\n", (*i)->Name.c_str(), (*i)->QueryBases.c_str());
    }
    fclose(foutput_mamf);
    if (is_after_aligned) {
        foutput_mamf=fopen(string(result_folder+"/"+gene_id+".mamf.single.reads").c_str(),"w+");
    }
    else {
        foutput_mamf=fopen(string(result_folder+"/"+gene_id+".mamf.before.single.reads").c_str(),"w+");
    }
    for (auto i = single_reads.begin(); i != single_reads.end(); i++) {
        fprintf(foutput_mamf, "%s\t%s\n", (*i)->Name.c_str(), (*i)->QueryBases.c_str());
    }
    fclose(foutput_mamf);
    FILE *finfo;
    finfo=fopen(string(result_folder +"/" + gene_id +".compatible.reads").c_str(), "w+");
    dump_compatible_reads(finfo);
    fclose(finfo);
}

void Earrings::test_memory_leak() {
    if (num_total_reads !=
        bam_reads.size() + unaligned_reads.size() + noused_reads.size())
        fprintf (stderr, "WARNING: Memory Leaking: total reads %zu\t record reads%zu\n", num_total_reads,  bam_reads.size() + unaligned_reads.size() + noused_reads.size());

}

int Earrings::load_read_data() {
    num_total_reads = 0;

    if (! jeweler_info->bam_reader.SetRegion(ref_id, left_pos, ref_id, right_pos)) {
        fprintf(stdout, "%s\n", jeweler_info->bam_reader.GetErrorString().c_str());
    }

    JewelerAlignment *al=new JewelerAlignment();

    while(jeweler_info->bam_reader.GetNextAlignment(*al)) {
        bam_reads.push_back(al);
        al=new JewelerAlignment();
    }
    delete al; // delete the last unused one
    num_total_reads = bam_reads.size();

    fprintf(stdout, "totally %zu bamalignments are loaded\n", num_total_reads);

    return 0;
}

int Earrings::load_transcript_data(bool is_prepare) {
    size_t i, j;

    transcript_helper(maternal_transcripts ,
                       jeweler_info->maternal_fasta, "maternal.",
                       result_folder + "/maternal.unmapped.bam",
                       is_prepare);

    transcript_helper(paternal_transcripts ,
                      jeweler_info->paternal_fasta, "paternal.",
                      result_folder + "/paternal.unmapped.bam",
                      is_prepare);
    if (is_prepare) {
        // fprintf(stdout, "%s\n", (string("python pipeline/mergeBam.py ")
        //                       + result_folder + "/maternal.unmapped.bam " +
        //                       result_folder + "/paternal.unmapped.bam "+
        //                       result_folder+" unmapped.bam").c_str());

        // system((string("python pipeline/mergeBam.py ")
        //      + result_folder + "/maternal.unmapped.bam " +
        //      result_folder + "/paternal.unmapped.bam " +
        //      result_folder+" unmapped.bam").c_str());
        system(string("cp " + result_folder + "/paternal.unmapped.bam " +
                      result_folder+"/unmapped.bam").c_str());

    }
    if (paternal_transcripts.size() != maternal_transcripts.size()
        || paternal_transcripts.size() == 0) {
        fprintf(stderr,
                "ERROR: number of transcripts does not match or no reads at all for gene %s  at %s:%d\n",
                gene_id.c_str(), __FILE__, __LINE__);
        exit(0);
    }

    for (i = 0; i < paternal_transcripts.size(); i++) {
        Transcript *p = paternal_transcripts[i];
        Transcript *m = maternal_transcripts[i];
        p->transcript_file = "paternal." + p->transcript_id + ".fasta";
        m->transcript_file = "maternal." + m->transcript_id + ".fasta";

        if (p->seq.size() != m->seq.size()) {
            fprintf(stderr,
                    "ERROR: transcript sequence size does not match at gene %s  at %s:%d\n",
                    p->transcript_id.c_str(),__FILE__, __LINE__);
            exit(0);
        }

        // find SNP position by given both paternal and maternal
        // transcripts sequences
        vector<unsigned int> snp_pos;
        vector<char> alleles;
        for (j = 0; j < p->seq.size(); j++) {
            if (p->seq[j] != m->seq[j]) {
                snp_pos.push_back(j);
                alleles.push_back(p->seq[j]);
            }
        }

        p->load_snps(snp_pos, alleles, EXON_PATERNAL);
        m->load_snps(snp_pos, alleles, EXON_MATERNAL);

        // mismatcher initialization
        mismatcher->add_transcript(p, TRANSCRIPT_PATERNAL);
        mismatcher->add_transcript(m, TRANSCRIPT_MATERNAL);
    }
    mismatcher->initialize();
    return 0;
}

template<class T>
vector<T *> duplicate_vector(vector<T *> in) {
    vector<T *> ret(in.size());
    size_t i;
    for (i = 0; i < in.size(); i ++) {
        ret[i] = new T(*in[i]);
    }
    return ret;
}

int Earrings::transcript_helper(vector<Transcript *> &transcripts,
                                FastaReference *fasta_ref,
                                string prefix,
                                string output_bam,
                                bool is_prepare) {
    //assuming gtf file has the same order of transcripts with the seq files
    size_t j;

    transcripts = duplicate_vector(jeweler_info->gene_id2transcripts[gene_id]);
    // TODO: put these code into gtf.cpp files.
    // All transcript class operations should be done in transcripts.
    TranscriptomeAligner ta;
    vector<string> reference_sequences;
    vector<string> reference_ids;
    string merged_fasta_file= (result_folder+"/"+prefix + transcripts[0]->gene_id+".fasta");
    for (j = 0; j < transcripts.size(); j++) {
        transcripts[j]->load_seq(fasta_ref);
        string filename = string(prefix+transcripts[j]->transcript_id);
        reference_sequences.push_back(result_folder + "/"+ filename);
        reference_ids.push_back(transcripts[j]->transcript_id);
    }
    if (is_prepare)
        ta.align(jeweler_info->left_unmapped_file, jeweler_info->right_unmapped_file,
                 reference_sequences, reference_ids, merged_fasta_file,output_bam);

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

            if (maternal_transcripts[j]->is_compatible(bam_reads[i], Transcript::tolerate)) {
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
void Earrings::dump_compatible_reads(FILE * fd) {
    size_t i;
    for (i = 0; i < compatible_reads.size(); i ++) {
        if (multiple_reads.find(compatible_reads[i]) ==
            multiple_reads.end())
            continue;
        fprintf(fd, "%s\t%zu", compatible_reads[i]->Name.c_str(),
                compatible_reads[i]->genome_position.size());

        for (size_t j = 0; j < compatible_reads[i]->genome_position.size(); j ++) {
            fprintf(fd, "\t%d", compatible_reads[i]->genome_position[j]);
        }
        fprintf(fd, "\n");
    }
}


void Earrings::align_reads() {
    size_t i,j,k;

	int num_maternal_alleles;
	int num_paternal_alleles;
	int total_alleles;
	bool is_compatible=false;
	set<JewelerAlignment *> cleared;
	vector<set<JewelerAlignment *> >noninfo;
	vector<int> paternal_alleles;
	vector<int> maternal_alleles;

	vector<JewelerAlignment *> new_bam_reads;
	vector<set<JewelerAlignment *> > read_lists;

	AlignmentGlue ag;
	Transcript::tolerate = 0;

	ag.glue(bam_reads, new_bam_reads, noused_reads);
	printf("%zu new paired reads, %zu unused reads\n", new_bam_reads.size(), noused_reads.size());
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

				maternal_transcripts[j]->register_read(bam_reads[i]);
				paternal_transcripts[j]->register_read(bam_reads[i]);

				maternal_transcripts[j]->match_alleles(bam_reads[i],
													   total_alleles,
													   &maternal_matcher);

				paternal_transcripts[j]->match_alleles(bam_reads[i],
													   total_alleles,
													   &paternal_matcher);

				//fprintf(stdout,"%d\n",total_alleles);
				num_maternal_alleles=maternal_alleles.size();
				num_paternal_alleles=paternal_alleles.size();
				if (total_alleles>0) {
					if (num_paternal_alleles == num_maternal_alleles) {
						noninfo[j].insert(bam_reads[i]);
						mismatcher->add_mismatches(maternal_transcripts[j],
												   bam_reads[i],
												   &maternal_matcher);
					}
					else {
						cleared.insert(bam_reads[i]);
						if (num_maternal_alleles > num_paternal_alleles) {
							maternal_transcripts[j]->register_allele_read(bam_reads[i]);
							mismatcher->add_mismatches(maternal_transcripts[j],
													   bam_reads[i],
													   &maternal_matcher);
						}
						else{
							paternal_transcripts[j]->register_allele_read(bam_reads[i]);
							mismatcher->add_mismatches(paternal_transcripts[j],
													   bam_reads[i],
													   &paternal_matcher);

						}
					}
				}
				else {
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
	FILE * finfo;

	fprintf(stdout, "%d compatible reads in %zu reads\n",
			num_compatible_reads,
			bam_reads.size()
			);

	if (sm!=NULL) {
		//count_multiple_alignments(/* is after alignment*/ true);
	}
	// finfo=fopen(string(result_folder+"/"+gene_id+".landscape.plot.meta").c_str(),"w+");

	// for (i=0;i<maternal_transcripts.size();i++) {
	// 	FILE *foutput;
	// 	foutput=fopen(string(result_folder+"/"+maternal_transcripts[i]->transcript_id+".landscape.plot.info").c_str(),"w+");
	// 	if (foutput == NULL) {
	// 		fprintf(stderr, "cannot open file at %s\n",
	// 				string(result_folder+"/"+maternal_transcripts[i]->transcript_id+".landscape.plot.info").c_str());
	// 	}
	// 	string name=maternal_transcripts[i]->transcript_id;
	// 	PileupPlot lp(maternal_transcripts[i],paternal_transcripts[i],noninfo[i], multiple_reads);
	// 	lp.generate_pileup_plot(finfo, foutput);
	// 	fclose(foutput);
	// }
	// fclose(finfo);

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
    classify_reads();
    for (auto i = multiple_reads.begin(); i != multiple_reads.end(); i++) {
        (*i)->dump_data(ed->add_multiple_read());
    }
    for (auto i = single_reads.begin(); i != single_reads.end(); i++) {
        (*i)->dump_data(ed->add_single_read());
    }
    return ;
}
