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

            num = self.load_num_reads(gene_id)
            if (num == 0 ) :
                print(gene_id)
                self.blacklist.add(gene_id)
        return result
        
    def __init__(self):
        cuffcompare_folder = sys.argv[1]
        self.jeweler_folder = sys.argv[2]
        self.bracelet_folder = sys.argv[3]
        self.mismatch_analyzer_folder = sys.argv[4]
        result_folder = sys.argv[5]
        print(os.path.dirname(result_folder))
        self.blacklist = set()
        cuffcompare_file = cuffcompare_folder + "/" + "cuffcompare.tracking"
        cuffcompare_result = self.load_cuffcompare_data(cuffcompare_file)

        fd = open(result_folder + "/blacklist.obj", "wb")
        pickle.dump(self.blacklist, fd)
        
    
if __name__ == "__main__":
    sg = SharedGraph()
