import re
import networkx as nx
import os
import json
import sys
import traceback
import argparse
import pydot
import math
from subprocess import Popen
from subprocess import call
import pickle
import constants
import shop_info
import cuffcompare
import jeweler_pb2
import zleveldb

class SharedGraph():
    num_genes = 0
    gene2gene = dict()

    def add_node_color(self, graph, name):
        if (self.is_unknown(name)):
            graph.add_node(name, color="black")
        elif (self.is_pseudo(name)):
            graph.add_node(name, color="red")
        else :
            graph.add_node(name, color="blue")

    def is_gene(self, name):
        return not self.is_pseudo(name) and not self.is_unknown(name)

    def is_unknown(self, name):
        cd = self.cuffcompare_result.get(name)
        assert cd, str(cd)
        return cd.is_unknown

    def is_pseudo(self, name):
        cd = self.cuffcompare_result.get(name)
        assert cd, str(cd)
        return cd.is_pseudo

    def register_edge(self, first, second):
        ## every edge is added twice
        d = (first, second)
        if self.is_gene(first):
            if self.is_gene(second):
                self.gene2gene[d] = constants.GG
            if self.is_pseudo(second):
                self.gene2gene[d] = constants.GP
            if self.is_unknown(second):
                self.gene2gene[d] = constants.GU
        elif self.is_pseudo(first):
            if self.is_gene(second):
                self.gene2gene[d] = constants.PG
            if self.is_pseudo(second):
                self.gene2gene[d] = constants.PP
            if self.is_unknown(second):
                self.gene2gene[d] = constants.PU
        elif self.is_unknown(first):
            if self.is_gene(second):
                self.gene2gene[d] = constants.UG
            if self.is_pseudo(second):
                self.gene2gene[d] = constants.UP
            if self.is_unknown(second):
                self.gene2gene[d] = constants.UU

    def load_mismatch_analyzer(self):
        data = jeweler_pb2.TranscriptMismatcherAnalyzerData()
        zleveldb.load_single_protobuf_data(self.shop_info.mismatch_analyzer_file,
                                           data)
        ret = dict()
        for d in data.location_result:
            if (d.coverage <7):
                continue
            if (d.num_mismatches < 5):
                continue
            if (d.pvalue > 0 and math.log10(d.pvalue)>-19):
                continue
            if (d.gene_id not in ret):
                ret[d.gene_id] = set()
            ret[d.gene_id].add(d.genome_location)
        return ret

    def load_bracelet_data(self, bracelet_file):
        graphs = list()
        graph_index = dict()
        data = jeweler_pb2.BraceletData()
        fd = open(bracelet_file, "rb")
        while zleveldb.load_protobuf_data(fd, data):
            print data.name
            print len(data.related_transcript)
            self.num_genes += 1
            graph = graph_index.get(data.name, None)
            if data.num_read == 0:
                self.black_list.add(data.name)
            if not graph:
                graph = nx.Graph()
                graphs.append(graph)
            for gene in data.related_transcript:
                graph.add_edge(data.name, gene.name,
                               label = gene.num_shared_read)
                print gene.num_shared_read
                self.add_node_color(graph, gene.name)
                self.add_node_color(graph, data.name)
                self.register_edge(data.name,
                                   gene.name)

    def __init__(self):
        self.shop_info = shop_info.ShopInfo()
        self.cuffcompare_result = \
            cuffcompare.CuffcompareResult(self.shop_info.cuffcompare_file)
        self.mismatch_analyzer = \
            self.load_mismatch_analyzer()
        self.bracelet_result = \
            self.load_bracelet_data(self.shop_info.bracelet_file)

        print("there are totally " +
              str(len(self.bracelet_result)) +
              " graphs are built .")
        print("there are totally " + str((self.num_genes)) +  " genes")
        pickle.dump(self.blacklist, open(self.blacklist_file, "wb"))
        # dump_dot_graph()

    def dump_dot_graph():
        fd = open(self.result_folder +"/" + "shared_graph","w+")
        k = 1
        for graph in self.bracelet_result:
            if (len(graph.nodes())>0):
                folder =  result_folder + "/graph." +str(k)
                dotfile = folder + "/graph.dot"
                psfile = folder + "/graph.png"
                if (not os.path.exists(folder)):
                    os.system("mkdir "+ folder)
                nx.write_dot(graph, dotfile)
                os.system("dot -Tpng " + dotfile  + " -o " + psfile)
                fd.write(folder+ "\t")
                for node in graph.nodes():
                    fd.write(node+"\t")
                fd.write("\n")
                k = k + 1

if __name__ == "__main__":
    sg = SharedGraph()
