#include "aligner.hpp"
 TranscriptomeAligner::TranscriptomeAligner(){
}

void TranscriptomeAligner::run_command(string command){
	fprintf(stdout, "%s  2> stderr_log \n", command.c_str());
	system((command+" 2> stderr_log\n").c_str());
}

void TranscriptomeAligner::align(string &left, string &right, vector<string>  &reference_sequences, vector<string> &reference_id, string & merged_fasta_file, string &output_file){
	string left_unmapped_file = left;
	string right_unmapped_file = right;
	ref_file = reference_sequences;
	vector<string> output_bams;
	string last, new_last;
	string command;
	run_command(string("samtools faidx " + merged_fasta_file));
	for (size_t i = 0; i < reference_sequences.size(); i ++){
		string left_sai = reference_sequences[i] + ".left.sai";
		string right_sai = reference_sequences[i] + ".right.sai";
		string output_bam = reference_sequences[i] + ".bam";
		string output_sam = reference_sequences[i] + ".sam";
		
		command = ((string("bwa index ") + reference_sequences[i]));
		run_command(command);
		command = ((string("bwa aln -n 10  ") + reference_sequences[i] 
					+ " " + left  + " > " + 
					left_sai).c_str());
		run_command(command);
		command = ((string("bwa aln -n 10 ") + reference_sequences[i] + " " + right  + " > " + 
					right_sai).c_str());
		run_command(command);

		command = ((string("bwa sampe ") + reference_sequences[i] + " " + 
					left_sai + " " + right_sai + " " + left + " " + right +  " 2> stderr_log | grep -E \"[1-9][0-9]*M\" > " + output_sam).c_str());
		run_command(command.c_str());		
		run_command(string("rm ") + left_sai);
		run_command(string("rm ") + right_sai);
		run_command(string("samtools view -bt " + merged_fasta_file + ".fai " + output_sam + " > " + output_bam));
		output_bams.push_back(output_bam);
	}
	if (output_bams.size() ==0) return ;
	command = "lib/bamtools/bin/bamtools merge ";
	for (size_t i=0; i < output_bams.size(); i ++){
		command += " -in " + output_bams[i];
	}
	command += " -out " + output_file;
	run_command(command);

}
