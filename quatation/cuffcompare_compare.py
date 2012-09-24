import re
import networkx as nx
import os
import json
import sys
import traceback
import argparse
import  pydot
import filecmp
import pickle
from subprocess import Popen
from subprocess import call


def load_appraiser_stat(app_file):
    lines = open(app_file).readlines()
    for line in lines:
        data = (line.strip().split('\t'))
        if (data[0] =="Total_alignments"):
            return data[1]

    return "0"

def load_cuffcompare_stat(cuffcompare_file, sample_id, filtered_gene_file, afd):
    lines = open( cuffcompare_file).readlines()
    result = dict()
    filtered_genes = set()
    if (filtered_gene_file is not None and os.path.exists(filtered_gene_file)):
        filtered_genes = pickle.load(open(filtered_gene_file, "rb"))
    transcripts = set()
    genes = set()
    genes_precision = set()
    unknown_transcripts = set()
    pseudo_transcripts = set()
    nonpseudo_transcripts = set()
    exact_transcripts = set()
    correct = 0
    correct_result_set = set()
    correct_result_set_precision = set()
    correct_transcript_result_set = set()
    correct_set = get_correct_data("../data/simulation_bam_new_merge/"+sample_id[0:9]+".abundance")
    correct_transcript_set = get_correct_transcript_data("../data/simulation_bam_new_merge/"+sample_id[0:9]+".abundance")

    pseudo_set = get_pseudogene_list()
    
    for line in lines:
        data = (line.strip().split('\t'))

        gene_id = (data[4][3:(data[4].find('|'))])

        transcript_id = data[4].split('|')[1]
        eid = ""
        tid = ""
        if (data[2].find("|") > 0):
            eid = data[2].split("|")[0]

        if (data[2].find("|") > 0):
            tid = data[2].split("|")[1]

        if (gene_id in filtered_genes):
            continue
        if (correct_set is not None and eid in correct_set):
            correct += 1
            correct_result_set.add(eid)
            correct_result_set_precision.add(gene_id)
        if (correct_transcript_set is not None and tid in correct_transcript_set):
            correct += 1
            correct_transcript_result_set.add(tid)
        
        transcripts.add(transcript_id)
        genes_precision.add(gene_id)
        if (eid != ""):
            genes.add(eid)
        else:
            genes.add(gene_id)
        if (data[2] != '-'):
            if (eid not in pseudo_set):
                nonpseudo_transcripts.add(eid)
                if (data[3] == '=' or data[3] == 'j'  ):
                    exact_transcripts.add(data[2])
            else:
                pseudo_transcripts.add(data[2])
        else:
            unknown_transcripts.add(gene_id)
    print((correct))
    print(len(correct_result_set))
    print(len(correct_set))
    print(len(correct_transcript_result_set))
    print(len(correct_transcript_set))
    if (len(correct_transcript_set) != 0):
        precision = len(correct_result_set_precision) / len(genes_precision)
        recall = len(correct_result_set) / len(genes)
        f1 = 2 *precision *recall / (precision + recall)
        print(len(correct_result_set_precision))
        print(len(genes_precision))
        print(str(precision) + "\t" +  str(recall)  +"\t"+ str(f1), file = afd)
        print(str(precision) + "\t" +  str(recall)  +"\t"+ str(f1))
    num_transcripts = len(transcripts)
    num_unknown_transcripts = len(unknown_transcripts)
    num_pseudo_transcripts = len(pseudo_transcripts)
    num_nonpseudo_transcripts = len(nonpseudo_transcripts)
    num_exact_transcripts = len(exact_transcripts)

    return dict(num_exact_transcripts = num_exact_transcripts,
                num_nonpseudo_transcripts = num_nonpseudo_transcripts,
                num_unknown_transcripts =num_unknown_transcripts,
                num_transcripts = num_transcripts,
                num_pseudo_transcripts = num_pseudo_transcripts
    )

def get_pseudogene_list():
    ret =set ()
    lines = open("../info/pseudo_gene_names").readlines()
    for line in lines:
        data = line.strip()
        ret.add(data)

    return ret
    
def get_pseudogene_transcript_list():
    ret =set ()
    lines = open("../info/pseudo_transcript_ids").readlines()
    for line in lines:
        data = line.strip()
        ret.add(data)

    return ret
    
def get_correct_transcript_data(filename):
    print(filename)
    if (filename is None ):
        return set()

    if (not os.path.isfile(filename)):
        return set()

    ret = set()
    lines = open( filename).readlines()
    for line in lines:
        data = line.strip().split(" ")
        ret.add(data[0])

    return ret
def get_correct_data(filename):
    print(filename)
    if (filename is None ):
        return set()

    if (not os.path.isfile(filename)):
        return set()

    ret = set()
    lines = open( filename).readlines()
    for line in lines:
        data = line.strip().split(" ")
        ret.add(data[3])

    return ret
        
def main():


    cuff1_root = "../result/merged_list/cuffcompare/"
    cuff2_root = "../result/new_merged_list/new_cuffcompare/"
    cuff3_root = "../result/new_merged_list/new_cuffcompare/"

    sg2_root = "../result/new_merged_list/new_shared_graph/"
    app1_root = "../result/merged_list/appraiser/"
    app2_root = "../result/new_merged_list/new_appraiser/"
    cuff1_root = "../result/simulation/cuffcompare/"
    cuff2_root = "../result/new_simulation/new_cuffcompare/"
    sg2_root = "../result/new_simulation/new_shared_graph/"
    app1_root = "../result/simulation/appraiser/"
    app2_root = "../result/new_simulation/new_appraiser/"
    cmp = filecmp.dircmp(cuff1_root, cuff2_root)
    fd = open("cuff.stats","w+")
    afd1 = open("accurate.stats1","w+")
    afd2 = open("accurate.stats2","w+")
    print("name\tbefore.exact\tbefore.pseudo\tbefore.nonpseudo\tbefore.unknown\tbefore.total\tafter.exact\tafter.pseudo\tafter.nonpseudo\tafter.unknown\tafter.total\tratio.exact\tratio.pseudo\tratio.nonpseudo\tratio.unknown\tratio.total", file = fd)
    fda = open("app.stats","w+")
    print("name\tbefore.alignments\tafter.alignments\tratio.alignments", file = fda)
    test_set = {"FF_0001_F_merged"}
    for i in cmp.common_dirs:
    ##for i in test_set:
        cuff1_result = load_cuffcompare_stat ( cuff1_root + i + "/cuffcompare.tracking",
                                               i,
                                               None, afd1)
        print(cuff1_result)                                               
        cuff2_result = load_cuffcompare_stat ( cuff2_root + i + "/cuffcompare.tracking",
                                               i,
                                               sg2_root + i + "/filtered_gene.obj", afd2)
        print(cuff2_result)
        r=[0]*5
        j = 0

        for v, w in zip(cuff1_result.values(), cuff2_result.values()):
            r[j] = (w-v)/v
            j += 1
        print(i+"\t"+"\t".join(["%s" % e for e in cuff1_result.values()]) +"\t"
              +("\t".join(["%s" % e for e in cuff2_result.values()])) + "\t"
              +("\t".join(["%.4lf" % e for e in r])) 
              ,file = fd)
        ##app1_result = load_appraiser_stat(app1_root + i +"/log")
        ##app2_result = load_appraiser_stat(app2_root + i +"/log")
        ##print(i + "\t" + app1_result + "\t" + app2_result + "\t"+str((int(app2_result)-int(app1_result))/int(app1_result)), file = fda)
    

if __name__ == "__main__":
    main()
