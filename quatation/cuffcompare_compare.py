import re
import networkx as nx
import os
import json
import sys
import traceback
import argparse
import  pydot
import filecmp
from subprocess import Popen
from subprocess import call


def load_appraiser_stat(app_file):
    lines = open(app_file).readlines()
    for line in lines:
        data = (line.strip().split('\t'))
        if (data[0] =="Total_alignments"):
            return data[1]

    return "0"

def load_cuffcompare_stat(cuffcompare_file):
    lines = open( cuffcompare_file).readlines()
    result = dict()
    num_transcripts = 0
    num_unknown_transcripts = 0
    num_pseudo_transcripts = 0
    num_nonpseudo_transcripts = 0
    num_exact_transcripts = 0
    for line in lines:
        data = (line.strip().split('\t'))
        num_transcripts += 1
        gene_id = (data[4][3:(data[4].find('|'))])
        if (data[2] != '-'):
            if (re.search('Gm[0-9]+',data[2]) is None):
                num_nonpseudo_transcripts +=1
                if (data[3] == '='):
                    num_exact_transcripts += 1
            else:
                num_pseudo_transcripts += 1
        else:
            num_unknown_transcripts += 1

    return dict(num_exact_transcripts = num_exact_transcripts,
                num_nonpseudo_transcripts = num_nonpseudo_transcripts,
                num_unknown_transcripts =num_unknown_transcripts,
                num_transcripts = num_transcripts,
                num_pseudo_transcripts = num_pseudo_transcripts
    )


def main():
    cuff1_root = "../result/merged_list/cuffcompare/"
    cuff2_root = "../result/new_merged_list/cuffcompare/"
    app1_root = "../result/merged_list/appraiser/"
    app2_root = "../result/new_merged_list/appraiser/"
    cmp = filecmp.dircmp(cuff1_root, cuff2_root)
    fd = open("cuff.stats","w+")
    print("name\tbefore.exact\tbefore.pseudo\tbefore.nonpseudo\tbefore.unknown\tbefore.total\tafter.exact\tafter.pseudo\tafter.nonpseudo\tafter.unknown\tafter.total\tratio.exact\tratio.pseudo\tratio.nonpseudo\tratio.unknown\tratio.total", file = fd)
    fda = open("app.stats","w+")
    print("name\tbefore.alignments\tafter.alignments\tratio.alignments", file = fda)

    for i in cmp.common_dirs:
        cuff1_result = load_cuffcompare_stat ( cuff1_root + i + "/cuffcompare.tracking" )
        cuff2_result = load_cuffcompare_stat ( cuff2_root + i + "/cuffcompare.tracking" )

        r=[0]*5
        j = 0
        print(cuff1_result)
        print(cuff1_result.values())
        for v, w in zip(cuff1_result.values(), cuff2_result.values()):
            r[j] = (w-v)/w
            j += 1
        print(r)
        print(i+"\t"+"\t".join(["%s" % e for e in cuff1_result.values()]) +"\t"
              +("\t".join(["%s" % e for e in cuff2_result.values()])) + "\t"
              +("\t".join(["%.4lf" % e for e in r])) 
              ,file = fd)            
        app1_result = load_appraiser_stat(app1_root + i +"/log")
        app2_result = load_appraiser_stat(app2_root + i +"/log")
        print(i + "\t" + app1_result + "\t" + app2_result + "\t"+str((int(app2_result)-int(app1_result))/int(app1_result)), file = fda)
    

if __name__ == "__main__":
    main()