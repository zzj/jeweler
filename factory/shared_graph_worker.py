import os
import json
import sys
import traceback
import argparse
from subprocess import call

def get_cuffcompare_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+args.cuffcompare_folder
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_cufflinks_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+args.cufflinks_folder
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_jeweler_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+ args.jeweler_folder
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_bracelet_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+ args.bracelet_folder
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_shared_graph_folder(args, filename):
    resultfolder='result/'+os.path.basename(args.filelist)+ args.shared_graph_folder
    
    resultsubfolder=resultfolder+'/'+os.path.basename(filename.strip().replace('.bam',''))
    return  resultsubfolder

def get_mismatch_analyzer(args, filename):
    alias = os.path.basename(filename.strip().replace('.bam',''))
    if (args.is_new_cufflinks):
        ma_result_folder = 'result/'+os.path.basename(args.filelist)+'/new_mismatch_analyzer/' + alias + '/'
    else:
        ma_result_folder = 'result/'+os.path.basename(args.filelist)+'/mismatch_analyzer/' + alias + '/'                
    return ma_result_folder
    
def shared_graph_worker(args):
    files=open(args.filelist).readlines()
    for f in files:
        cuffcompare_folder=get_cuffcompare_folder(args, f)
        jeweler_folder=get_jeweler_folder(args, f)
        bracelet_folder=get_bracelet_folder(args,f)
        mismatch_analyzer_folder = get_mismatch_analyzer(args, f)
        shared_graph_folder = get_shared_graph_folder(args, f)
        sample_id =  os.path.basename(f.strip().replace('.bam',''))

        ##stupid python evoke a R program that cannot read a file
        print('python3.2 shop/shared_graph.py '+cuffcompare_folder + "/ " +jeweler_folder + "/ "+bracelet_folder + "/ " + mismatch_analyzer_folder + "/ "+
              shared_graph_folder + "/ " +sample_id  )

def classify_gene_worker(args):
    files=open(args.filelist).readlines()
    for f in files:
        shared_graph_folder = get_shared_graph_folder(args, f)
        cufflinks_folder = get_cufflinks_folder(args, f)
        sample_id =  os.path.basename(f.strip().replace('.bam',''))

        ##stupid python evoke a R program that cannot read a file
        print('python3.2 shop/classify_genes.py ' + shared_graph_folder + "/ " +sample_id +" " + cufflinks_folder )
        ##print ("~/bin/bin/R CMD BATCH --no-save --no-restore '--args name=\""+sample_id+"\"' shop/shared_graph.R temp/"+sample_id)
