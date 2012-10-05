"""
This object is converetd a cuffcompare file to a protobuf project and used for further analysis.
"""

import sys
import jeweler_pb2
import os
import zleveldb
import pseudo
import re

class CuffcompareResult:

    @property
    def protobuf_file(self):
        return self.filename + ".proto.data"

    def __init__(self, filename):
        self.filename = filename
        self.data = jeweler_pb2.CuffcompareData()
        self.pseudo_set = pseudo.PseudoSet()
        if False:
            zleveldb.load_single_protobuf_data(self.protobuf_file, self.data)
        else:
            self.load_cuffcompare_raw_data(self.data)
            zleveldb.write_single_protobuf_data(self.protobuf_file, self.data)

        self.index = dict()
        self.gene2num_exon = dict()
        self.transcript2num_exon = dict()
        self.load_gtf()
        self.build_index()

    def load_gtf(self):
        lines = open("data/database/Mus_musculus.NCBIM37.63.gtf").readlines()
        for line in lines:
            line_data = line.strip().split('\t')
            if line_data[2] != "exon":
                continue
            match = re.search(r'gene_name "(.*)"; t', line_data[8])
            gene_name = match.group(1)
            match = re.search(r'transcript_id "([0-9A-Z]*)', line_data[8])
            transcript_name = match.group(1)
            match = re.search(r'exon_number "([0-9]*)', line_data[8])
            num_exon = int(match.group(1))

    def load_cuffcompare_raw_data(self, data):
        lines = open(self.filename).readlines()

        for line in lines:
            gene_data = data.gene_data.add()
            line_data = line.strip().split('\t')

            gene_data.gene_id = line_data[4][3:(line_data[4].find('|'))]
            gene_data.transcript_id = line_data[4].split('|')[1]
            if line_data[2] != '-':
                gene_data.transcript_name = line_data[2].split("|")[1]
                gene_data.gene_name = line_data[2].split("|")[0]
            else:
                gene_data.is_unknown = True
            gene_data.fpkm = float(line_data[4].split("|")[6])
            gene_data.matched_score = line_data[3]
            if line_data[2] != '-' and (line_data[3] == '=' or line_data[3] == 'j'):
                gene_data.is_good = True

            if self.pseudo_set.is_pseudo_gene(gene_data.gene_name):
                 gene_data.is_pseudo = True
                 ## this is super important, should be unit test
                 gene_data.is_good = False

    def build_index(self):
        for d in self.data.gene_data:
            self.index[d.gene_id] = d

    def count(self):
        self.num_unknown_transcripts = 0
        self.num_pseudo_transcripts = 0
        self.num_nonpseudo_transcripts = 0
        self.num_good_transcripts = 0
        self.num_transcripts = len(self.data.gene_data)
        good_genes = set()
        for d in self.data.gene_data:
            if d.is_unknown:
                self.num_unknown_transcripts += 1
                continue
            if d.is_pseudo:
                self.num_pseudo_transcripts += 1
            else:
                self.num_nonpseudo_transcripts += 1
            if d.is_good:
                self.num_good_transcripts += 1

    def get_meta_count(self):
        self.count()
        return dict(num_good_transcripts = self.num_good_transcripts,
                    num_nonpseudo_transcripts = self.num_nonpseudo_transcripts,
                    num_unknown_transcripts = self.num_unknown_transcripts,
                    num_transcripts = self.num_transcripts,
                    num_pseudo_transcripts = self.num_pseudo_transcripts
        )


    def all_gene_names(self, black_list = None, only_good = False):
        gene_names = []
        for d in self.data.gene_data:
            if black_list and d.gene_id in black_list:
                continue
            gene_names.append(d.gene_name)
        return gene_names

    def all_transcript_names(self, black_list = None, only_good = False):
        transcript_names = []
        for d in self.data.gene_data:
            if black_list and d.gene_id in black_list:
                continue
            transcript_names.append(d.transcript_name)
        return transcript_names

    def get(self, gene_id):
        return self.index.get(gene_id, None)
