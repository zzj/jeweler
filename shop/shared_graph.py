import re
import networkx as nx
import os
import json
import sys
import traceback
import argparse
import  pydot
from subprocess import Popen
from subprocess import call
from gene_meta import GeneMeta

class SharedGraph:
    genes = dict()
    num_genes = 0
    num_pseudos = 0
    num_unknowns = 0
    
    gene2gene = dict()
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
    genelist = set()
    good_gene = set()
    def load_cuffcompare_data(self, cuffcompare_file):
        lines = open( cuffcompare_file).readlines()
        result = dict()
        for line in lines:
            data = (line.strip().split('\t'))
            gene_id = (data[4][3:(data[4].find('|'))])
            if (data[2] != '-'):
                result[gene_id] = gene_id + "|"+data[2]
            else:
                result[gene_id] = gene_id
            if (data[3] == '='):
                self.good_gene.add(gene_id + "|" + data[2])
        return result


    def add_node_color(self, graph, name):
        if (self.is_unknown(name)):
            graph.add_node(name, color="black")
        elif (self.is_pseudo(name)):
            graph.add_node(name, color="red")
        else :
            graph.add_node(name, color="blue")

    def is_pseudo(self, name):
        return (re.search('Gm[0-9]+',name) is not None)

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
        self.generate_training_data(d, self.gene2gene[d])
                                    
    def register_node(self, name, genename):
        if (name in self.genelist):
            return
        self.genelist.add(name)
        if (self.is_gene(name)):
            self.genes[name] = GeneMeta(genename, self.jeweler_folder, self.GENE)
            self.num_genes += 1
        if (self.is_pseudo(name)):
            self.genes[name] = GeneMeta(genename, self.jeweler_folder, self.PSEUDO)
            self.num_pseudos += 1
        if (self.is_unknown(name)):
            self.genes[name] = GeneMeta(genename, self.jeweler_folder, self.UNKNOWN)
            self.num_unknowns += 1

    def load_bracelat_data(self, bracelet_file, cuffcompare):
        lines = open( bracelet_file).readlines()
        result = dict()
        graphs = list()
        idx = 0
        print(len(lines))
        for line in lines:
            data = (line.strip().split('\t'))
            if (len(data)==2):
                continue
            graph = None
            idx += 1

            for i in range(int(len(data)/2)):
                if (data[i*2]  in result):
                    graph = result[data[i*2]]
            if (data[0] not in cuffcompare):
                cuffcompare[data[0]] = data[0]
                
            self.register_node(cuffcompare[data[0]], data[0])
            if (graph is None):
                graph=nx.Graph()
                graphs.append(graph)
            for i in range(int(len(data)/2-1)):
                if (int(data[i*2+3]) > 10):
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
        trainset = set([self.GG, self.GP, self.PG, self.PP])
        if (category not in trainset):
            return

        origin = relationship [0]
        target = relationship [1]
        if(category == self.GG):
            if (origin not in self.good_gene or target not in self.good_gene):
                return 
        if(category == self.GP):
            if (origin not in self.good_gene):
                return 
        if(category == self.PG):
            if ( target not in self.good_gene):
                return 
        filename = (self.bracelet_folder + "/" + self.genes[origin].genename
                    + "/" + self.genes[target].genename)
        lines = open(filename).readlines()
        filename = (self.bracelet_folder + "/" + self.genes[origin].genename
                    + "/" + self.genes[target].genename)
        locations1 = []
        coverage1 = []
        for line in lines:
            data = line.strip().split('\t')
            locations1.append(data[0])
            coverage1.append(int(data[1]))
        filename = (self.bracelet_folder + "/" + self.genes[target].genename
                    + "/" + self.genes[origin].genename)
        locations2 = []
        coverage2 = []
        lines = open(filename).readlines()
        for line in lines:
            data = line.strip().split('\t')
            locations2.append(data[0])
            coverage2.append(int(data[1]))
        print(str(category)+"\t" + 
              str(self.genes[origin].get_error_rate(locations1)) + "\t" +
              str(self.genes[target].get_error_rate(locations2)) + "\t" +
              str(self.genes[origin].get_shared_coverage(coverage1)) + "\t" +
              str(self.genes[target].get_shared_coverage(coverage2)) + "\t" +
              str(self.genes[origin].get_covered_ratio(locations1)) + "\t" +
              str(self.genes[target].get_covered_ratio(locations2)), file = self.fd)

        
    def __init__(self):
        cuffcompare_folder = sys.argv[1]
        self.jeweler_folder = sys.argv[2]
        self.bracelet_folder = sys.argv[3]
        result_folder = sys.argv[4]
        sample_id = sys.argv[5]
        print(os.path.dirname(result_folder))
        if (not os.path.exists(result_folder)):
            os.makedirs(result_folder)
        cuffcompare_file = cuffcompare_folder + "/" + "cuffcompare.tracking"
        cuffcompare_result = self.load_cuffcompare_data(cuffcompare_file)

        bracelet_file = self.bracelet_folder + "/" + "result.bracelet.1"
        print(bracelet_file)
        self.fd = open("training","w+")
        bracelet_result = self.load_bracelat_data(bracelet_file, cuffcompare_result)

        print("there are totally " + str(len(bracelet_result))+ " graphs are built .")
        print("there are totally " + str((self.num_genes)) +  " genes")
        print("there are totally " + str((self.num_pseudos)) +  " pseudos")
        print("there are totally " + str((self.num_unknowns)) +  " unknowns")
        ##self.generate_training_data()
        return 
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
