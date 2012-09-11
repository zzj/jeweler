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

class GeneMeta:
    coverage = dict()
    mismatches = dict()
    total_coverage = 0
    num_isoforms = 0
    num_exons = 0

    def __init__(self, genename, root, category, etype, pattern):
        self.genename = genename
        self.root = root
        self.category = category
        self.pattern = pattern
        self.expressing_type = etype
        self.load_data()

    def load_data(self):
        mismatcher_file = self.root+"/"+self.genename+"/"+self.genename+".mismatcher"
        lines = open(mismatcher_file).readlines()
        lineid = 1
        colname = dict()
        locationidx = 0
        coverageidx = 1
        mismatchesidx = 2
        for line in lines:
            data = (line.strip().split('\t'))
            if (lineid ==1):
                for i in range(0,len(data)):
                    colname[data[i]] = i

                if ('location' not in colname):
                    print("No location column\n", file = sys.stderr)
                    return 
                if ('coverage' not in colname):
                    print("No coverage column\n", file = sys.stderr)
                    return 
                if ('mismatches' not in colname):
                    print("No mismatches column\n", file = sys.stderr)
                    return
                coverageidx = colname['coverage']
                locationidx = colname['location']
                mismatchesidx = colname['mismatches']
            else:
                self.coverage[data[locationidx]] = int(data[coverageidx])
                self.mismatches[data[locationidx]] = int(data[mismatchesidx])
                self.total_coverage += int(data[coverageidx])
            lineid += 1
        print(self.genename + "\t" + str(len(self.coverage)))
        filename = (self.root 
                    + "/" + self.genename
                    + "/" + self.genename + ".landscape.plot.meta")
        lines = open(filename).readlines()
        self.num_isoforms = len(lines);
        num_exons1 = 0
        for line in lines:
            data = line.strip().split('\t')
            num_exons1 += int(data[1])
        self.num_exons = num_exons1 / self.num_isoforms
        

    def get_error_rate(self, locations):
        num_mismatches = 0
        total = 0.001
        for l in locations:
            if (l in self.coverage):
                num_mismatches += self.mismatches [l]
                total +=self.coverage[l]
        return num_mismatches/total

    def get_covered_ratio(self, locations):
        return float(len(locations))/(len(self.coverage)+0.01)

    def get_shared_coverage(self, cs):
        k = 0
        for c in cs:
            k += c
        return k/float(self.total_coverage +0.01)
