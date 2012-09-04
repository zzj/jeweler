import re
import networkx as nx
import os
import json
import sys
import traceback
import argparse
import  pydot
import math
from subprocess import Popen
from subprocess import call
from gene_meta import GeneMeta
from gene_relationship import GeneRelationship
import pickle 
import constants

class Classifier:

    def count_num_not_sure(self):
        ret = 0
        for gene, genemeta in self.genes.items():
            if (genemeta.expressing_type == constants.NOT_SURE):
                ret += 1
        return ret
                
    def __init__(self):
        self.root = sys.argv[1]
        self.sample_id = sys.argv[2]
        self.cufflinks_path = sys.argv[3]
        fd = open(sys.argv[1] + "/" + "gene_relationships.obj", "rb")
        self.gene_relationships = pickle.load(fd)
        fd = open(sys.argv[1] + "/" + "gene_meta.obj", "rb")
        self.genes = pickle.load( fd)
        fd = open(sys.argv[1] + "/" + "blacklist.obj", "rb")
        self.blacklist = pickle.load( fd)
        self.correct_set = self.get_correct_data("data/simulation_bam_new_merge/"+self.sample_id[0:9]+".abundance")
        self.num_not_sure = self.count_num_not_sure()
        print(str(self.num_not_sure))
        print(len(self.genes))
        print("finished loading")
        self.classify()
        
    def get_correct_data(self, filename):
        print(filename)
        if (filename is None ):
            return None

        if (not os.path.isfile(filename)):
            return None

        ret = set()
        lines = open( filename).readlines()
        for line in lines:
            data = line.strip().split(" ")
            ret.add(data[3])

        return ret

    def gene_rule_expressed_exact_match_gene(self, gene_meta):
        if ((gene_meta.pattern == "=" or gene_meta.pattern == "j") 
            and gene_meta.category == constants.GENE):
        ##if( gene_meta.category == self.GENE):
            return constants.EXPRESSED
        return constants.NOT_SURE
        
        
    def relationship_rule_unexpressed_cannot_share_cover_too_much(self, gene_relationship):
        if (gene_relationship.origin_covered_ratio <
            gene_relationship.target_covered_ratio - 0.3):
            return True
        return False
        
    def relationship_rule_unexpressed_cannot_share_read_too_much(self, gene_relationship):
        if (gene_relationship.origin_shared_coverage < 
            gene_relationship.target_shared_coverage - 0.2):
            return True
        return False
        
    def relationship_rule_unexpressed_cannot_beat_a_good_gene(self, gene_relationship):
        if (self.genes[gene_relationship.origin].num_exons > 1 and 
            self.genes[gene_relationship.target].num_exons == 1):
            return True
        return False

    def relationship_rule_unexpressed_has_bad_match(self, gene_relationship):
        if(gene_relationship.origin_consistent_locations < gene_relationship.target_consistent_locations):
            return True
        return False

    def get_gene_name(self, name):
        return self.genes[name].genename

    def get_gene_expressing_type(self, name):
        return self.genes[name].expressing_type

    def get_gene_category(self, name):
        return self.genes[name].category

    def count_gene_types(self, gene_set):
        num_genes = 0
        num_pseudos = 0
        num_unknowns = 0
        correct_set = set()
        wrong_set = set()
        for gene in  gene_set:
            genename = gene.strip().split("|")
            if self.get_gene_category(gene) != constants.GENE and \
               self.get_gene_category(gene) != constants.PSEUDO and \
               self.get_gene_category(gene) != constants.UNKNOWN:
                continue
                
            if (self.get_gene_category(gene) == constants.GENE):
                num_genes += 1
            if (self.get_gene_category(gene) == constants.PSEUDO):
                num_pseudos += 1
            if (self.get_gene_category(gene) == constants.UNKNOWN):
                num_unknowns += 1

            if (len(genename)==3 and self.correct_set is not None and genename[1] in self.correct_set):
                correct_set.add(gene)
            else:
                wrong_set.add(gene)

        print(str(num_genes), str(num_pseudos) ,str(num_unknowns), (len(correct_set)), (len(wrong_set)))

    def classify_normal_relationship(self, rule, filter_set ):
        ret = set()
        print(rule.__name__)
        for gene_pair, gene_relationship in self.gene_relationships.items():
            origin = gene_relationship.origin
            target = gene_relationship.target
            if ( self.get_gene_expressing_type(origin)!= constants.EXPRESSED or
                 self.get_gene_expressing_type(target)!= constants.NOT_SURE):
                continue
            if (target in filter_set):
                continue

            if (rule(gene_relationship)):
                ret.add(target)
        return ret
        
    def classify_good_relationship(self, rule, filter_set ):
        ret = set()
        print(rule.__name__)
        for gene_pair, gene_relationship in self.gene_relationships.items():
            origin = gene_relationship.origin
            target = gene_relationship.target
            if ( self.get_gene_expressing_type(origin)== constants.EXPRESSED or
                 self.get_gene_expressing_type(target)!= constants.NOT_SURE):
                continue
            if (target in filter_set):
                continue

            if (rule(gene_relationship)):
                ret.add(target)
        return ret
        
    def count_expressed_related_genes(self):
        study_set = set()
        gene_set = set()
        self.count_gene_types(self.genes.keys())
        gene_set = self.classify_good_relationship(self.relationship_rule_unexpressed_cannot_beat_a_good_gene, set())
        self.count_gene_types(gene_set)
        gene_set2 = self.classify_normal_relationship(self.relationship_rule_unexpressed_cannot_share_read_too_much, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        gene_set2 = self.classify_good_relationship(self.relationship_rule_unexpressed_has_bad_match, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        gene_set2 = self.classify_normal_relationship(self.relationship_rule_unexpressed_has_bad_match, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        gene_set2 = self.classify_good_relationship(self.relationship_rule_unexpressed_cannot_share_read_too_much, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        gene_set2 = self.classify_normal_relationship(self.relationship_rule_unexpressed_cannot_beat_a_good_gene, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        gene_set2 = self.classify_normal_relationship(self.relationship_rule_unexpressed_cannot_share_cover_too_much, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        gene_set2 = self.classify_good_relationship(self.relationship_rule_unexpressed_cannot_share_cover_too_much, gene_set)
        self.count_gene_types(gene_set2)
        gene_set = gene_set.union(gene_set2)
        self.count_gene_types(gene_set)
        gene_names = set()

        for r in gene_set:
            gene_names.add(self.get_gene_name(r))
        gene_names = gene_names.union(self.blacklist)
        pickle.dump(gene_names, open(self.root+"/filtered_gene.obj", "wb"))

    def dump_gtf_file(self):
        gene_names = pickle.load(open(self.root+"/filtered_gene.obj", "rb"))
        lines = open(self.cufflinks_path+"/"+"transcripts.gtf").readlines()
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
        open(self.cufflinks_path+"/"+"new_transcripts.gtf","w+").writelines(new_lines)
        new_lines = filter(is_filtered, lines)
        open(self.cufflinks_path+"/"+"suspicious_transcripts.gtf","w+").writelines(new_lines)
        transcripts_lines = [line for line in lines if is_filtered(line) and is_transcript(line)]
        fd = open(self.cufflinks_path+"/"+"sus_distribution.txt","w+")

        for line in transcripts_lines:
            k = line.split('\t')
            chrid = (k[0][3:])
            start = (k[3])
            end = (k[4])
            print(chrid + "\t"+ start +  "\t"+ end , file = fd)
        fd.close()

    def classify(self):
        num_expressed = 0
        for gene, genemeta in self.genes.items():
            genemeta.expressing_type = self.gene_rule_expressed_exact_match_gene(genemeta)
            if (genemeta.expressing_type == constants.EXPRESSED):
                num_expressed += 1
        print("Found out " + str(num_expressed) + " expressed genes")
        self.count_expressed_related_genes()
        self.dump_gtf_file()


if __name__ == "__main__":
    
    classifier = Classifier()
