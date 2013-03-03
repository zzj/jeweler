from __future__ import print_function
import sys
sys.path.insert(0, 'shop/')

import re
import networkx as nx
import os
import json


import traceback
import argparse
import pydot
import filecmp
import pickle
from subprocess import Popen
from subprocess import call
import cuffcompare


def load_cuffcompare_stat(cuffcompare_file):
    return cuffcompare.CuffcompareResult(cuffcompare_file).get_meta_count()

def main():
    datafolder = "merged_list"
    datafolder = "inbreds"
    cuff1_root = "result/" + datafolder + "/cuffcompare/"
    cuff2_root = "result/new_" + datafolder + "/new_cuffcompare/"
    cmp = filecmp.dircmp(cuff1_root, cuff2_root)
    fd = open("quatation/" + datafolder + ".cuff.stats","w+")
    print("name\tbefore.unknown\tbefore.pseudo\tbefore.nonpseudo\tbefore.exact\tbefore.total\tafter.unknown\tafter.pseudo\tafter.nonpseudo\tafter.exact\tafter.total\tratio.unknown\tratio.pseudo\tratio.nonpseudo\tratio.exact\tratio.total", file = fd)
    for i in cmp.common_dirs:
        cuff1_file = cuff1_root + i + "/cuffcompare.tracking"
        cuff2_file = cuff2_root + i + "/new_cuffcompare.tracking"
        if not (os.path.exists(cuff1_file) and os.path.exists(cuff2_file)):
            continue
        cuff1_result = load_cuffcompare_stat(cuff1_file)
        print(cuff1_result)
        cuff2_result = load_cuffcompare_stat(cuff2_file)
        print(cuff2_result)
        r= [0] * 5
        j = 0

        for v, w in zip(cuff1_result.values(), cuff2_result.values()):
            r[j] = (w-v) * 1.0/v
            j += 1
        print(i+"\t"+"\t".join(["%s" % e for e in cuff1_result.values()]) +"\t"
              +("\t".join(["%s" % e for e in cuff2_result.values()])) + "\t"
              +("\t".join(["%.4lf" % e for e in r]))
              ,file = fd)

if __name__ == "__main__":
    main()
