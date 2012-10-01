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
import cuffcompare
import pprint

from sklearn import svm
from sklearn import tree
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from sklearn.metrics import zero_one_score
from sklearn.neighbors.nearest_centroid import NearestCentroid
from sklearn.ensemble import RandomForestClassifier

def f1(precision, recall):
    return 2 * precision * recall / (precision + recall)


class JewelerClassifier:
    def __init__(self):
        self.shop_info = shop_info.ShopInfo()
        self.cuffcompare_result = \
            cuffcompare.CuffcompareResult(self.shop_info.cuffcompare_file)
        self.training_data = pickle.load(open(self.shop_info.test_data_file, "rb"))
        self.black_list = pickle.load(open(self.shop_info.blacklist_file, "rb"))
        self.bad_genes = set()
        self.good_genes = set()
        self.bad_transcripts = set()
        self.good_transcripts = set()
        self.fit = None
        print(self.shop_info.sample_id[0:9])
        self.correct_gene_set, self.correct_transcript_set = self.get_correct_data("data/simulation_bam_new_merge/" + self.shop_info.sample_id[0:9] + ".abundance")
        self.generate_training_data()
        self.train()
        self.classify()
        self.dump_gtf_file()

    def get_correct_data(self, filename):
        if filename is None :
            return (None, None)

        if not os.path.isfile(filename):
            return (None, None)

        gene = set()
        transcript =set()
        lines = open(filename).readlines()
        for line in lines:
            data = line.strip().split(" ")
            gene.add(data[3])
            transcript.add(data[0])
        return (gene, transcript)

    def train(self):
        raise NotImplementedError

    def print_output(self, goods, correct_set):
        pre = len(goods.intersection(correct_set)) * 1.0 / len(goods)
        recall = len(goods.intersection(correct_set)) * 1.0 / len(correct_set)
        print (pre, recall, f1(pre, recall))

    def classify(self):
        predict = self.fit.predict(self.X)
        print(sum(self.Y) * 1.0 / len(self.Y))
        print(zero_one_score(self.Y, predict))
        print(precision_score(self.Y, predict))
        print(recall_score(self.Y, predict))
        print(f1_score(self.Y, predict))

        for i, v in enumerate(predict):
            if v == 1:
                self.black_list.add(self.training_data[self.idx[i]].target_name)
                # if self.training_data[self.idx[i]].target_gene_name in self.correct_gene_set:
                    # print (self.training_data[self.idx[i]].origin_gene_name,
                    #        self.training_data[self.idx[i]].target_gene_name)
                    # print (self.training_data[self.idx[i]].X)
        if self.correct_gene_set:
            all_genes = self.cuffcompare_result.all_gene_names()
            self.good_genes = self.cuffcompare_result.all_gene_names(self.black_list)
            self.print_output(self.good_genes, self.correct_gene_set)
            self.print_output(all_genes, self.correct_gene_set)

            all_transcripts = self.cuffcompare_result.all_transcript_names()
            self.good_transcripts = self.cuffcompare_result.all_transcript_names(self.black_list)
            self.print_output(self.good_transcripts, self.correct_transcript_set)
            self.print_output(all_transcripts, self.correct_transcript_set)

        predict = self.fit.predict(self.allX)

        for i, v in enumerate(predict):
            if v == 1:
                self.black_list.add(self.training_data[i].target_name)

    def dump_gtf_file(self):
        gene_names = self.black_list
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
        self.allX = []
        self.allY = []
        self.idx = []
        for i, s in enumerate(self.training_data):
            if s.origin_name in self.black_list or s.target_name in self.black_list:
               continue
            self.allX.append(s.X)
            self.allY.append(0)
            if s.used_for_training(self.correct_gene_set, self.black_list):
                self.X.append(s.X)
                self.Y.append(s.Y(self.correct_gene_set))
                self.idx.append(i)

        pickle.dump((self.X, self.Y),
                    open("real_data/" + self.shop_info.sample_id, "wb"))

class SVMJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = svm.SVC(class_weight = {0:1, 1:1})
        ## http://scikit-learn.org/0.11/auto_examples/plot_precision_recall.html
        half = len(self.X)
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])
        pickle.dump(self.fit, open("learning_model", "wb"))

class TreeJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = tree.DecisionTreeClassifier()
        ## http://scikit-learn.org/0.11/auto_examples/plot_precision_recall.html
        half = len(self.X)
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])
        pickle.dump(self.fit, open("learning_model", "wb"))

class RFJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = RandomForestClassifier()
        ## http://scikit-learn.org/0.11/auto_examples/plot_precision_recall.html
        half = len(self.X)
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])
        pickle.dump(self.fit, open("learning_model", "wb"))

class KNNJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = NearestCentroid()
        ## http://scikit-learn.org/0.11/auto_examples/plot_precision_recall.html
        half = len(self.X) / 2
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])

class ClassifyJewelerClassifier(JewelerClassifier):
    def train(self):
        self.fit = pickle.load(open("learning_model", "rb"))

if __name__ == "__main__":
    sjc = ClassifyJewelerClassifier()
