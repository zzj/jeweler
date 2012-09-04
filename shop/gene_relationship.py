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

class GeneRelationship:
    def __init__(self, shared_graph, relationship, category, mismatch_analyzer):
        self.category = category
        self.origin = relationship[0]
        origin = relationship[0]
        self.target = relationship[1]
        target = relationship[1]
        self.is_origin_good = False
        self.is_target_good = False
        self.origin_consistent_locations = 0
        self.target_consistent_locations = 0
        if (origin in shared_graph.good_gene):
            self.is_origin_good = True
        if ( target in shared_graph.good_gene):
            self.is_target_good = True
        origin_genename = shared_graph.genes[origin].genename
        target_genename = shared_graph.genes[target].genename
        filename = (shared_graph.bracelet_folder
                    + "/" + origin_genename
                    + "/" + target_genename)
        lines = open(filename).readlines()
        locations1 = []
        coverage1 = []
        consistent_locations = set()
        if (origin_genename in mismatch_analyzer):
            consistent_locations = mismatch_analyzer[origin_genename]

        for line in lines:
            data = line.strip().split('\t')
            locations1.append(int(data[0]))
            coverage1.append(int(data[1]))
            if (int(data[0]) in consistent_locations):
                self.origin_consistent_locations += 1
        filename = (shared_graph.bracelet_folder
                    + "/" + target_genename
                    + "/" + origin_genename)
        locations2 = []
        coverage2 = []
        lines = open(filename).readlines()
        consistent_locations = set()
        if (target_genename in mismatch_analyzer):
            consistent_locations = mismatch_analyzer[target_genename]
        for line in lines:
            data = line.strip().split('\t')
            locations2.append(int(data[0]))
            coverage2.append(int(data[1]))
            if (int(data[0]) in consistent_locations):
                self.target_consistent_locations += 1
        
        self.overlap = (len(set(locations1) & set(locations2)))
        
        if (shared_graph.genes[origin].get_covered_ratio(coverage1)> 1):
            print("Error: " + origin + " does not correct")
            print(len(shared_graph.genes[origin].coverage))
            print(len(coverage1))
        
        self.origin_error_rate = shared_graph.genes[origin].get_error_rate(locations1)
        self.target_error_rate = shared_graph.genes[target].get_error_rate(locations1)
        self.origin_shared_coverage = shared_graph.genes[origin].get_shared_coverage(coverage1)
        self.target_shared_coverage = shared_graph.genes[target].get_shared_coverage(coverage1)
        self.origin_covered_ratio = shared_graph.genes[origin].get_covered_ratio(coverage1)
        self.target_covered_ratio = shared_graph.genes[target].get_covered_ratio(coverage1)
        self.origin_num_isoforms = shared_graph.genes[origin].num_isoforms
        self.target_num_isoforms = shared_graph.genes[target].num_isoforms
        self.origin_num_exons = shared_graph.genes[origin].num_exons
        self.target_num_exons = shared_graph.genes[target].num_exons

        
        
        
        