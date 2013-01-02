#include <iostream>
#include <string>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamConstants.h>
#include <Fasta.h>
#include "laboratory/cigar_holder.hpp"
#include "main/jeweler_alignment.hpp"
using namespace BamTools;
using namespace std;


int main (int argc, char *argv[]) {
	BamReader reader;
	JewelerAlignment al;

	FastaReference fr;


	if (!reader.Open(argv[1])) {
		fprintf(stderr,"Cannot open bam file!\n");
		exit(0);
	}

	fr.open(argv[2]);
	auto references=reader.GetReferenceData();

	for (size_t i=0;i<references.size();i++) {
		fprintf(stdout,"%s\n",references[i].RefName.c_str());
	}

	fprintf(stdout,"fasta\n");
	for (size_t i= 0;
		 i<fr.index->sequenceNames.size();
		 i++) {
		fprintf(stdout,"%s\n",fr.index->sequenceNames[i].c_str());
	}
	while (reader.GetNextAlignment(al)) {

		string refseq=fr.getSubSequence(references[al.RefID].RefName,
										al.Position,
										al.Length);
		if (!al.IsReverseStrand()) continue;
		fprintf(stdout,"%s\n",refseq.c_str());
		fprintf(stdout,"%s\n",al.QueryBases.c_str());
		fprintf(stdout,"%s\n",get_cigar_string(al).c_str());
		fprintf(stdout,"%s\n",al.TagData.c_str());
		uint16_t nm;
		al.GetTag("NM",nm);
		fprintf(stdout,"%d\n",nm);


	}

	return 0;
}
