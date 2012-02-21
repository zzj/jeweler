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


def load_cuffcompare_data(cuffcompare_file):
    lines = open( cuffcompare_file).readlines()
    result = dict()
    for line in lines:
        data = (line.strip().split('\t'))
        gene_id = (data[4][3:(data[4].find('|'))])
        if (data[2] != '-'):
            result[gene_id] = gene_id + "|"+data[2]
        else:
            result[gene_id] = gene_id

    return result


def add_node_color(graph, name):
    if (name.find("|")<0):
        graph.add_node(name, color="black")
    elif (re.search('Gm[0-9]+',name) is None):
        graph.add_node(name, color="red")
    else :
        graph.add_node(name, color="blue")

def load_bracelat_data(bracelet_file, cuffcompare):
    lines = open( bracelet_file).readlines()
    result = dict()
    graphs = list()
    print(len(lines))
    for line in lines:
        data = (line.strip().split('\t'))
        if (len(data)==2):
            continue
        graph = None
        
        for i in range(int(len(data)/2)):
            if (data[i*2]  in result):
                graph = result[data[i*2]]
                
        if (graph is None):
            graph=nx.Graph()
            graphs.append(graph)
        for i in range(int(len(data)/2-1)):
            if (int(data[i*2+3]) > 10):
                result[data[i*2+2]] = graph
                if (data[0] not in cuffcompare):
                    cuffcompare[data[0]] = data[0]
                if (data[i*2+2] not in cuffcompare):
                    cuffcompare[data[i*2+2]] = data[i*2+2]
                graph.add_edge(cuffcompare[data[0]],
                               cuffcompare[data[i*2+2]],weight = int(data[i*2+3]))
                add_node_color(graph, cuffcompare[data[0]])
                add_node_color(graph, cuffcompare[data[i*2+2]])
                
    return graphs

def main():
    cuffcompare_folder = sys.argv[1]
    jeweler_folder = sys.argv[2]
    bracelet_folder = sys.argv[3]
    result_folder = sys.argv[4]
    sample_id = sys.argv[5]
    print(os.path.dirname(result_folder))
    if (not os.path.exists(result_folder)):
        os.makedirs(result_folder)
    cuffcompare_file = cuffcompare_folder + "/" + "cuffcompare.tracking"
    cuffcompare_result = load_cuffcompare_data(cuffcompare_file)

    bracelet_file = bracelet_folder + "/" + "result.bracelet"
    bracelet_result = load_bracelat_data(bracelet_file, cuffcompare_result)
    
    print("there are totally " + str(len(bracelet_result))+ " graphs are built .")
    fd = open(result_folder +"/" + "shared_graph","w+")
    k=1


    
    for graph in bracelet_result:
        if (len(graph.nodes())>0):
            folder =  result_folder + "/graph." +str(k)
            dotfile = folder + "/graph.dot"
            psfile = folder + "/graph.ps"
            if (not os.path.exists(folder)):
                os.system("mkdir "+ folder)
            nx.write_dot(graph, dotfile)
            os.system("dot -Tps " + dotfile  + " -o " + psfile)
            fd.write(folder+ "\t")
            for node in graph.nodes():
                fd.write(node+"\t")
            fd.write("\n")
            k = k + 1

    
if __name__ == "__main__":
    main()
