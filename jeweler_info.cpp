#include "jeweler_info.hpp"
#include "transcript.hpp"
#include "constants.hpp"
#include "jeweler_alignment.hpp"

int JewelerInfo::check_args(const int i, char *argv[], const char * name, string &a) {
	if (strcmp(argv[i], name) == 0) {
		a = argv[i + 1];
		return i + 1;
	}
	return i;
}

vector<Transcript *> JewelerInfo::get_paternal_transcripts(const string &gene_id) {
    return this->gene_id2paternal_transcripts[gene_id];
}

vector<Transcript *> JewelerInfo::get_maternal_transcripts(const string &gene_id) {
    return this->gene_id2maternal_transcripts[gene_id];
}

int JewelerInfo::get_num_genes() {
    return this->gene_id2maternal_transcripts.size();
}


int JewelerInfo::build_gene_id2transcripts() {
	fprintf(stdout, "Bulding gene_id to transcripts index ...\n");
	gene_id2maternal_transcripts.clear();
	gene_id2paternal_transcripts.clear();
	for (size_t i = 0; i < maternal_transcripts.size(); i++) {
		gene_id2maternal_transcripts[maternal_transcripts[i]->gene_id()].push_back(maternal_transcripts[i]);
		gene_id2paternal_transcripts[paternal_transcripts[i]->gene_id()].push_back(paternal_transcripts[i]);
	}
	gene_id.clear();
	for (auto j = gene_id2maternal_transcripts.begin();
		 j != gene_id2maternal_transcripts.end();
		 j ++) {
		gene_id.push_back(j->first);
	}
	return 0;
}

int JewelerInfo::get_refID(string chr) {
	for (size_t i = 0; i < references.size(); i++) {
		if (references[i].RefName == chr) {
			return i;
		}
	}
	return NOT_FOUND;
}

JewelerInfo::JewelerInfo(int argc, char *argv []) {
	int i;
	for (i = 1; i < argc; i ++) {
		i = check_args(i, argv, "-maternal_strain_ref_file", maternal_strain_ref_file);
		i = check_args(i, argv, "-paternal_strain_ref_file", paternal_strain_ref_file);
		i = check_args(i, argv, "-maternal_strain_id", maternal_strain_id);
		i = check_args(i, argv, "-paternal_strain_id", paternal_strain_id);
		i = check_args(i, argv, "-bam_file", bam_file);
		i = check_args(i, argv, "-alias", alias);
		i = check_args(i, argv, "-result_file", result_file);
		i = check_args(i, argv, "-result_folder", result_folder);
		i = check_args(i, argv, "-left_unmapped_file", left_unmapped_file);
		i = check_args(i, argv, "-right_unmapped_file", right_unmapped_file);
		i = check_args(i, argv, "-gtf_input_file", gtf_input_file);
	}
	result_folder += string("/") + alias + "/";
	create_directory(path(result_folder));
	
	maternal_fasta = new FastaReference();
	paternal_fasta = new FastaReference();
	maternal_fasta->open(maternal_strain_ref_file);
	paternal_fasta->open(paternal_strain_ref_file);

	load_gtf_file(gtf_input_file, maternal_transcripts, paternal_transcripts,
                  maternal_fasta, paternal_fasta);
	this->build_gene_id2transcripts();
	fprintf(stdout, "There are totally %zu transcripts loaded from %s\n", 
			maternal_transcripts.size(), gtf_input_file.c_str());
    open_bam(this->bam_reader, this->bam_file);
	references = bam_reader.GetReferenceData();
}

JewelerInfo::JewelerInfo() {

}

JewelerInfo::~JewelerInfo() {

}

