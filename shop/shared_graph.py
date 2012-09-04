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
from gene_meta import GeneMeta
from gene_relationship import GeneRelationship
import pickle 

class SharedGraph:
    genes = dict()
    num_genes = 0
    num_suspicious_pseudos = 0
    num_suspicious_unknowns = 0
    num_unknowns = 0
    num_pseudos = 0

    gene2gene = dict()
    gene_relationships = dict()
    GG = "GG"
    GP = "GP"
    GU = "GU"
    PG = "PG"
    PP = "PP"
    PU = "PU"
    UG = "UG"
    UP = "UP"
    UU = "UU"
    GENE = 0
    PSEUDO = 1
    UNKNOWN = 2
    NOT_SURE = -1
    EXPRESSED = 0
    UNEXPRESSED = 1

    genelist = set()
    good_gene = set()
    gene_pattern = dict()
    suspicious_pseudo_genes = set()
    suspicious_unknown_genes = set()

    def load_num_reads(self, gene_id):
        ret = 0
        single = 0
        multiple = 0
        filename  = self.jeweler_folder + "/"+ gene_id + "/" + gene_id + ".mamf.before.single.reads"
        lines = open( filename).readlines()
        single += len(lines)
        filename  = self.jeweler_folder + "/"+ gene_id + "/" + gene_id + ".mamf.before.multiple.reads"
        lines = open( filename).readlines()
        multiple += len(lines)
        filename  = self.jeweler_folder + "/"+ gene_id + "/" + gene_id + ".mamf.multiple.reads"
        lines = open( filename).readlines()
        ret += len(lines)
        filename  = self.jeweler_folder + "/"+ gene_id +"/"+gene_id+ ".mamf.single.reads"
        lines = open( filename).readlines()
        ret += len(lines)
        return ret

    def load_cuffcompare_data(self, cuffcompare_file):
        lines = open( cuffcompare_file).readlines()
        result = dict()
        for line in lines:
            data = (line.strip().split('\t'))
            gene_id = (data[4][3:(data[4].find('|'))])
            tid = ""
            if (data[2] != '-'):
                result[gene_id] = gene_id + "|"+data[2]
                tid = data[2].split("|")[1]
            else:
                result[gene_id] = gene_id

            self.gene_pattern[gene_id] = data[3]
            if (data[3] == '='):
                self.good_gene.add(gene_id + "|" + data[2])

            if (self.is_pseudo(result[gene_id])):
                self.pseudos[tid] = 0

            if (self.is_unknown(result[gene_id])):
                self.unknowns.add(result[gene_id])
        return result

    def get_pseudogene_list(self):
        ret =set ()
        lines = open("info/pseudo_gene_names").readlines()
        for line in lines:
            data = line.strip()
            ret.add(data)
        return ret

    def add_node_color(self, graph, name):
        if (self.is_unknown(name)):
            graph.add_node(name, color="black")
        elif (self.is_pseudo(name)):
            graph.add_node(name, color="red")
        else :
            graph.add_node(name, color="blue")

    def is_pseudo(self, name):
        eid = "none"
        tid = "none"
        if (name.find("|") > 0):
            eid = name.split("|")[1]
            tid = name.split("|")[2]
        if (eid in self.pseudo_set):
            return True
        else:
            return False

    def is_unknown(self, name):
        return (name.find("|")<0)

    def is_gene(self, name):
        return ( not self.is_pseudo(name) and  not self.is_unknown(name))

    def register_edge(self, first, second):
        ## every edge is added twice
        d = (first, second)
        if (self.is_gene(first)):
            if (self.is_gene(second)):
                self.gene2gene[d] = self.GG
            if (self.is_pseudo(second)):
                self.gene2gene[d] = self.GP
            if (self.is_unknown(second)):
                self.gene2gene[d] = self.GU
        elif (self.is_pseudo(first)):
            if (self.is_gene(second)):
                self.gene2gene[d] = self.PG
            if (self.is_pseudo(second)):
                self.gene2gene[d] = self.PP
            if (self.is_unknown(second)):
                self.gene2gene[d] = self.PU
        elif (self.is_unknown(first)):
            if (self.is_gene(second)):
                self.gene2gene[d] = self.UG
            if (self.is_pseudo(second)):
                self.gene2gene[d] = self.UP
            if (self.is_unknown(second)):
                self.gene2gene[d] = self.UU
        self.gene_relationships[d] =  self.generate_training_data(d, self.gene2gene[d])
                                    
    def register_node(self, name, genename):
        if (name in self.genelist):
            return
        self.genelist.add(name)
        if (self.is_gene(name)):
            self.genes[name] = GeneMeta(genename, self.jeweler_folder, self.GENE,
                                        self.NOT_SURE, self.gene_pattern[genename])
            self.num_genes += 1
        if (self.is_pseudo(name)):
            self.pseudos[name.split("|")[2]] += 1
            self.genes[name] = GeneMeta(genename, self.jeweler_folder, self.PSEUDO,
                                        self.NOT_SURE, self.gene_pattern[genename])
            self.num_suspicious_pseudos += 1
        if (self.is_unknown(name)):
            self.genes[name] = GeneMeta(genename, self.jeweler_folder, self.UNKNOWN,
                                        self.NOT_SURE, self.gene_pattern[genename])
            self.num_suspicious_unknowns += 1

    def load_mismatch_analyzer(self, mismatch_analyzer_file):
        ret = dict()
        lines = open( mismatch_analyzer_file).readlines()
        id = -1
        for line in lines:
            data = (line.strip().split('\t'))
            id += 1
            if (id == 0):
                continue

            total = int(data[4])
            miss = int(data[5])
            p_value = float(data[2])
            gene_id = data[0]
            loc = int(data[1])
            if (total <7):
                continue
            if (miss < 5):
                continue
            if (p_value > 0 and math.log10(p_value)>-19):
                continue
            if (gene_id not in ret):
                ret[gene_id] = set()
                
            ret[gene_id].add(loc)
            
        return ret
        

    def load_bracelat_data(self, bracelet_file, cuffcompare):
        lines = open( bracelet_file).readlines()
        result = dict()
        graphs = list()
        idx = 0
        for line in lines:
            data = line.strip().split('\t')
            if len(data)==2:
                continue
            graph = None
            idx += 1

            for i in range(int(len(data)/2)):
                if data[i*2] in result:
                    graph = result[data[i*2]]
            if (data[0] not in cuffcompare):
                cuffcompare[data[0]] = data[0]
                
            self.register_node(cuffcompare[data[0]], data[0])
            if (graph is None):
                graph=nx.Graph()
                graphs.append(graph)
            for i in range(int(len(data)/2-1)):
                if (int(data[i*2+3]) > 0):
                    result[data[i*2+2]] = graph
                    if (data[i*2+2] not in cuffcompare):
                        cuffcompare[data[i*2+2]] = data[i*2+2]
                    graph.add_edge(cuffcompare[data[0]],
                                   cuffcompare[data[i*2+2]],label = int(data[i*2+3]))
                    self.register_node(cuffcompare[data[i*2+2]], data[i*2+2])
                    self.register_edge(cuffcompare[data[0]],
                                  cuffcompare[data[i*2+2]])
                    self.add_node_color(graph, cuffcompare[data[0]])
                    self.add_node_color(graph, cuffcompare[data[i*2+2]])
        return graphs

    def generate_training_data(self, relationship, category):
        return GeneRelationship(self, relationship, category, self.mismatch_analyzer)
        
    def __init__(self):
        cuffcompare_folder = sys.argv[1]
        self.jeweler_folder = sys.argv[2]
        self.bracelet_folder = sys.argv[3]
        self.mismatch_analyzer_folder = sys.argv[4]
        result_folder = sys.argv[5]
        print(os.path.dirname(result_folder))
        self.pseudo_set = self.get_pseudogene_list()
        self.pseudos = dict()
        self.unknowns = set()
        self.blacklist = set()
        if (not os.path.exists(result_folder)):
            os.makedirs(result_folder)
        cuffcompare_file = cuffcompare_folder + "/" + "cuffcompare.tracking"
        cuffcompare_result = self.load_cuffcompare_data(cuffcompare_file)
        mismatch_analyzer_file = self.mismatch_analyzer_folder + "/" + "result.consistent.locations"
        self.mismatch_analyzer = self.load_mismatch_analyzer(mismatch_analyzer_file)

        bracelet_file = self.bracelet_folder + "/" + "result.bracelet"

        print(bracelet_file)
        self.fd = open("training","w+")
        bracelet_result = self.load_bracelat_data(bracelet_file, cuffcompare_result)
        for k, v in self.pseudos.items():
            if (v == 0):
                print(k)

        self.num_pseudos = len(self.pseudos)
        self.num_unknowns = len(self.unknowns)
        print(len(self.blacklist))
        print("there are totally " + str(len(bracelet_result))+ " graphs are built .")
        print("there are totally " + str((self.num_genes)) +  " genes")
        print("there are totally " +
              str((self.num_suspicious_pseudos)) + "/" + str(self.num_pseudos)+
              " pseudos")
        print("there are totally " +
              str((self.num_suspicious_unknowns)) + "/" + str(self.num_unknowns) +
              " unknowns")
        ##self.generate_training_data()
        fd = open(result_folder + "/gene_relationships.obj", "wb")
        pickle.dump(self.gene_relationships, fd)
        fd = open(result_folder + "/gene_meta.obj", "wb")
        pickle.dump(self.genes, fd)

        fd = open(result_folder +"/" + "shared_graph","w+")
        k=1

        for graph in bracelet_result:
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
