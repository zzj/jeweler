from __future__ import print_function
import re
import networkx as nx
import os
import json
import sys
import traceback
import argparse
import pydot
import math
from gene_relationship import GeneRelationship
import pickle
import constants
import shop_info

from sklearn import svm
from sklearn import tree
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

class JewelerClassifier:
    def __init__(self):
        self.shop_info = shop_info.ShopInfo()
        self.training_data = pickle.load(open(self.shop_info.test_data_file, "rb"))
        self.black_list = pickle.load(open(self.shop_info.blacklist_file, "rb"))
        self.bad_genes = set()
        self.good_genes = set()
        self.generate_training_data()
        self.classify()
        self.dump_gtf_file()

    def classify(self):
        raise NotImplementedError

    def dump_gtf_file(self):
        gene_names = self.bad_genes
        lines = open(self.shop_info.cufflinks_folder + "/" +
                     "transcripts.gtf").readlines()

        def is_filtered(line):
            start = line.find("gene_id \"")
            end = line.find("\";", start)
            return line[(start+9):end] in gene_names

        def is_transcript(line):
            start = line.find("exon")
            return start < 0

        def is_not_filtered(line):
            return not is_filtered(line)

        new_lines = filter(is_not_filtered, lines)
        open(self.shop_info.cufflinks_folder+"/"+"new_transcripts.gtf","w+").writelines(new_lines)
        new_lines = filter(is_filtered, lines)
        open(self.shop_info.cufflinks_folder+"/"+"suspicious_transcripts.gtf","w+").writelines(new_lines)
        transcripts_lines = [line for line in lines if is_filtered(line) and is_transcript(line)]
        fd = open(self.shop_info.cufflinks_folder+"/"+"sus_distribution.txt","w+")
        for line in transcripts_lines:
            k = line.split('\t')
            chrid = (k[0][3:])
            start = (k[3])
            end = (k[4])
            print(chrid + "\t"+ start +  "\t"+ end , file = fd)
        fd.close()

    def generate_training_data(self):
        self.X = []
        self.Y = []
        for s in self.training_data:
            if s.origin_name in self.black_list or s.target_name in self.black_list:
               continue
            if s.used_for_training:
                self.X.append(s.X)
                self.Y.append(s.Y)


class SVMJewelerClassifier(JewelerClassifier):
    def classify(self):
        clf = svm.SVC(kernel='linear', probability=True)
        ## http://scikit-learn.org/0.11/auto_examples/plot_precision_recall.html
        half = len(self.X) / 2
        probas_ = clf.fit(self.X[0:half], self.Y[0:half]).predict_proba(self.X[half:])
        precision, recall, thresholds = precision_recall_curve(Y[half:], probas_[:, 1])
        area = auc(recall, precision)
        print("Area Under Curve: %0.2f" % area)
        print (zip(precision, recall, thresholds))


class TreeJewelerClassifier(JewelerClassifier):
    def classify(self):
        clf = tree.DecisionTreeClassifier()
        ## http://scikit-learn.org/0.11/auto_examples/plot_precision_recall.html
        half = len(self.X) / 2
        probas_ = clf.fit(self.X[0:half], self.Y[0:half]).predict_proba(self.X[half:])
        precision, recall, thresholds = precision_recall_curve(self.Y[half:],
                                                               probas_[:, 1])
        area = auc(recall, precision)
        print("Area Under Curve: %0.2f" % area)
        print (zip(precision, recall, thresholds))


if __name__ == "__main__":
    sjc = SVMJewelerClassifier()
