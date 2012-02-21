import os
import json
import sys
import traceback
import argparse
from subprocess import call

def get_cuffcompare_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+'/cuffcompare'
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_jeweler_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+'/jeweler'
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_bracelet_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+'/bracelet'
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_shared_graph_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+'/shared_graph/'
    
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def shared_graph_worker(args):
    files=open(args.filelist).readlines()
    for f in files:
        cuffcompare_folder=get_cuffcompare_folder(args, f)
        jeweler_folder=get_jeweler_folder(args, f)
        bracelet_folder=get_bracelet_folder(args,f)
        shared_graph_folder = get_shared_graph_folder(args, f)
        sample_id =  os.path.basename(f.strip().replace('.bam',''))

        ##stupid python evoke a R program that cannot read a file
        ##print('python3.2 shop/shared_graph.py '+cuffcompare_folder + "/ " +jeweler_folder + "/ "+bracelet_folder + "/ " + shared_graph_folder + "/ " +sample_id + )
        print ("~/bin/bin/R CMD BATCH --no-save --no-restore '--args name=\""+sample_id+"\"' shop/shared_graph.R temp/"+sample_id)
