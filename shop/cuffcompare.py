"""
This object is converetd a cuffcompare file to a protobuf project and used for further analysis.
"""

import sys
import jeweler_pb2
import os
import zleveldb
import pseudo

class CuffcompareResult:

    @property
    def protobuf_file(self):
        return self.filename + ".proto.data"

    def __init__(self, filename):
        self.filename = filename
        self.data = jeweler_pb2.CuffcompareData()
        self.pseudo_set = pseudo.PseudoSet()
        if os.path.exists(self.protobuf_file):
            zleveldb.load_single_protobuf_data(self.protobuf_file, self.data)
        else:
            self.load_cuffcompare_raw_data(self.data)
            zleveldb.write_single_protobuf_data(self.protobuf_file, self.data)

        self.index = dict()
        self.build_index()

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

            gene_data.matched_score = line_data[3]
            if line_data[2] != '-' and line_data[3] == '=':
                gene_data.is_good = True
                
            if self.pseudo_set.is_pseudo_gene(gene_data.gene_name):
                 gene_data.is_pseudo = True

    def build_index(self):
        for d in self.data.gene_data:
            self.index[d.gene_id] = d

    def get(self, gene_id):
        return self.index.get(gene_id, None)
