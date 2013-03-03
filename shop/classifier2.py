from __future__ import print_function
import re
import networkx as nx
import os
import json
import sys
import traceback
import collections
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
    if precision + recall < 0.000000001:
        return 0
    return 2 * precision * recall / (precision + recall)

class JewelerClassifier:
    def __init__(self):
        self.shop_info = shop_info.ShopInfo()
        self.cuffcompare_result = \
            cuffcompare.CuffcompareResult(self.shop_info.cuffcompare_file)
        self.training_data = \
            pickle.load(open(self.shop_info.test_data_file, "rb"))
        self.black_list = \
            pickle.load(open(self.shop_info.blacklist_file, "rb"))
        self.fit = None
        simulation_file = "data/simulation_bam_new_merge/" + \
                          self.shop_info.sample_id[0:9] + ".abundance"
        self.correct_gene_set, self.correct_transcript_set = \
             self.get_correct_data(simulation_file)
        self.generate_training_data()
        self.train()
        self.classify()
        self.dump_gtf_file()
        self.dump_cuffcompare_file()

    def get_correct_data(self, filename):
        if filename is None:
            return (None, None)

        if not os.path.isfile(filename):
            return (None, None)

        gene = []
        transcript = []
        lines = open(filename).readlines()
        for line in lines:
            data = line.strip().split(" ")
            gene.append(data[3])
            transcript.append(data[0])
        return (gene, transcript)

    def train(self):
        raise NotImplementedError


    ## the following three functions can be moved out from this class.
    def calculate_gene_accuracy(self, goods, alls, correct_set):
        '''Calculate the accuracy per gene.

        If the gene is occurred multiple times in the result, only count once as correct one.
        '''
        goods_counter = collections.Counter(goods)
        correct_counter = collections.Counter(correct_set)
        num_correct_genes = sum((goods_counter & correct_counter).values())
        pre = num_correct_genes * 1.0 / len(alls)
        correct_number = len(filter(set(goods).__contains__, set(correct_set))) * 1.0
        recall = correct_number / len(set(correct_set))
        return (pre, recall , f1(pre, recall))

    def calculate_error_rate(self, bads, alls, correct_set):
        '''Calculate the error rate per gene'''
        num_bad_genes = len([i for i in bads if i not in set(correct_set)])
        pre = float(num_bad_genes) / len(bads)
        residue = set(alls) - set(correct_set)
        recall = float(len(bads & residue)) / len(residue)
        total = float(len(bads)) / len(correct_set)
        return(pre, recall, total)

    def calculate_testing_accuracy(self, Y, predict):
        redundancy = sum(self.Y) * 1.0 / len(self.Y)
        accuracy = zero_one_score(Y, predict)
        precision = precision_score(Y, predict)
        recall = recall_score(Y, predict)
        f1 = f1_score(Y, predict)

    def classify(self):
        predict = self.fit.predict(self.X)
        self.calculate_testing_accuracy(self.Y, predict)
        for i, v in enumerate(predict):
            if v == 1:
                self.black_list.add(self.training_data[self.idx[i]].target_name)
                # if self.training_data[self.idx[i]].target_gene_name in self.correct_gene_set:
                    # print (self.training_data[self.idx[i]].origin_gene_name,
                    #        self.training_data[self.idx[i]].target_gene_name)
                    # print (self.training_data[self.idx[i]].X)
        predict = self.fit.predict(self.allX)

        for i, v in enumerate(predict):
            if v == 1:
                self.black_list.add(self.training_data[i].target_name)
        if self.correct_gene_set:
            self.old_cuffcompare_result = \
                  cuffcompare.CuffcompareResult("result/simulation/cuffcompare/" + self.shop_info.sample_id + "/cuffcompare.tracking")
            # self.old_cuffcompare_result = \
            #       cuffcompare.CuffcompareResult("result/simulation/cuffcompare/" + self.shop_info.sample_id + "/cuffcompare.tracking")
            all_genes = self.cuffcompare_result.all_gene_names()
            good_genes = self.cuffcompare_result.all_gene_names(self.black_list)
            old_all_genes = self.old_cuffcompare_result.all_gene_names()
            old_good_genes = self.old_cuffcompare_result.all_gene_names()
            bad_genes = set(all_genes) - set(good_genes)

            print(self.calculate_gene_accuracy(good_genes, all_genes,
                                               self.correct_gene_set))
            print(self.calculate_gene_accuracy(old_good_genes, old_all_genes,
                                               self.correct_gene_set))
            print(self.calculate_error_rate(bad_genes, all_genes,
                                            self.correct_gene_set))

            old_good_transcripts = \
                 self.old_cuffcompare_result.all_transcript_names(self.black_list)
            old_all_transcripts = self.old_cuffcompare_result.all_transcript_names()
            all_transcripts = \
                 self.cuffcompare_result.all_transcript_names(self.black_list)
            good_transcripts = \
                self.cuffcompare_result.all_transcript_names(self.black_list)
            self.calculate_gene_accuracy(good_transcripts, all_transcripts,
                              self.correct_transcript_set)
            self.calculate_gene_accuracy(old_good_transcripts, old_all_transcripts,
                              self.correct_transcript_set)

    def dump_cuffcompare_file(self):
        gene_names = self.black_list
        lines = open(self.shop_info.cuffcompare_file).readlines()
        def cuffcompare_filter(line):
            line_data = line.strip().split('\t')
            gene_id = line_data[4][3:(line_data[4].find('|'))]
            return gene_id not in gene_names
        new_lines = filter(cuffcompare_filter, lines)
        open(self.shop_info.cuffcompare_folder+"/"+"new_cuffcompare.tracking","w+").writelines(new_lines)

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
        transcripts_lines = [line for line in lines
                             if is_filtered(line) and is_transcript(line)]
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
                    open("simulation_data/" + self.shop_info.sample_id, "wb"))

class SVMJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = svm.SVC(class_weight = {0:1, 1:1})
        half = len(self.X)
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])
        pickle.dump(self.fit, open("svm_learning_model", "wb"))

class TreeJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = tree.DecisionTreeClassifier()
        half = len(self.X)
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])
        pickle.dump(self.fit, open("tree_learning_model", "wb"))

class RFJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = RandomForestClassifier()
        half = len(self.X)
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])
        pickle.dump(self.fit, open("rf_learning_model", "wb"))

class KNNJewelerClassifier(JewelerClassifier):
    def train(self):
        clf = NearestCentroid()
        half = len(self.X) / 2
        self.fit = clf.fit(self.X[0:half], self.Y[0:half])

class ClassifyJewelerClassifier(JewelerClassifier):
    def train(self):

        # remove the current instance from the training set.
        #data = "simulation_data/" + self.shop_info.sample_id[0:9] + "_merged"
        #os.system("python shop/generate_learning_model.py " + data)

        self.fit = pickle.load(open("learning_model_all", "rb"))

if __name__ == "__main__":
    sjc = ClassifyJewelerClassifier()
