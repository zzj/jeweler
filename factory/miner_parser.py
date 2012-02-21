
import os
import json
import sys
import traceback
import argparse
from subprocess import call

def initialize_parser():
        parser=argparse.ArgumentParser(description=
                                       'The main batch script for RNA-seq analysis.')
        
        parser.add_argument('--json',
                            help='legendary json function. Process one json file a time',
                            dest='json')

        parser.add_argument('--filelist',
                            help='A list of bam files, for multipe bamfiles. ',
                            dest='filelist')

        parser.add_argument('--reftable',
                            help='A reference table for reference genomes',
                            dest='reftable')

        parser.add_argument('--result_folder',
                            help='result folder',
                            dest='result_folder')

        parser.add_argument('--alias',
                            help='project alias',
                            dest='alias')

        parser.add_argument('--maternal_strain_ref_seq',
                            help='matran strain fasta file',
                            dest='maternal_strain_ref_seq')

        parser.add_argument('--maternal_strain_id',
                            help='matran strain id',
                            dest='maternal_strain_id')

        parser.add_argument('--paternal_strain_ref_seq',
                            help='matran strain fasta file',
                            dest='paternal_strain_ref_seq')

        parser.add_argument('--paternal_strain_id',
                            help='matran strain id',
                            dest='paternal_strain_id')

        parser.add_argument('--bam_file',
                            help='bam files',
                            dest='bam_file')

        parser.add_argument('--gtf_input_file',
                            help='gtf input file',
                            dest='gtf_input_file')

        parser.add_argument('--result_file',
                            help='result files',
                            dest='result_file');


        parser.add_argument('--cuffcompare',
                            action='store_true',
                            dest='is_cuffcompare')

        parser.add_argument('--cufflinks',
                            action='store_true',
                            dest='is_cufflinks')

        parser.add_argument('--jeweler',
                            action='store_true',
                            dest='is_jeweler')

        parser.add_argument('--appraiser',
                            action='store_true',
                            dest='is_appraiser')

        parser.add_argument('--jeweler_only',
                            action='store_true',
                            dest='is_jeweler_only')

        parser.add_argument('--mismatch_analyzer',
                            action='store_true',
                            dest='is_mismatch_analyzer')
        
        parser.add_argument('--bracelet',
                            action='store_true',
                            dest='is_bracelet')

        parser.add_argument('--shared_graph',
                            action='store_true',
                            dest='plot_shared_graph')
                            
        return parser
