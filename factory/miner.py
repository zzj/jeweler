import os
import json
import sys
import traceback
import argparse
from subprocess import call


def split_gff_files(left_strain_ref_seq, right_strain_ref_seq,
					left_strain_id, right_strain_id,
					bam_file,
					result_folder, output_info_file, result_file):
	result_file_fd= open(result_file,"w+")
	try:
		os.makedirs()
	except:
		print("Notice: the folder was created.\n")


	with open(configuration['gtf_input_file']) as f:
		gene_meta = f.readlines()

	last_id		  = ""
	chr			 = ""
	left_region	 = -1
	right_region	= -1
	is_file_created = False
	tranfile		= None
	id			  = 0

	for line in gene_meta:
		info=line.strip().split('\t')
		if (len(info)!=9):
			break
		# Find a new start of a transcript
		if (info[2]=='transcript'):
			gene_id=info[8].split(';')[0].split('"')[1]

			if (last_id==gene_id):
				# new transcript, but same gene with previous one
				chr=info[0]
				if ((left_region)<0):
					left_region = int(info[3])
				else :
					left_region = min(int(info[3]),left_region)

				if (right_region<0):
					right_region=int(info[4])
				else :
					right_region=max(int(info[4]),right_region)
			else :
				# new transcript, new gene
				if (tranfile!=None):
					tranfile.close()
					id=id+1
					print(str(id)+"\n")
					left_output_seq  = gene_folder+last_id+'.'+left_strain_id+".fa"
					right_output_seq = gene_folder+last_id+'.'+right_strain_id+".fa"
					read_seq		 = gene_folder+last_id+".seq.bam";
					call(['gffread', '-w', left_output_seq, '-g', left_strain_ref_seq, 
						  gene_folder+last_id])
					call(['gffread', '-w', right_output_seq, '-g', right_strain_ref_seq, 
						  gene_folder+last_id])
					call(["samtools","view", bam_file , 
						 chr +":" +str(left_region)  +"-" +str(right_region) ,
						 "-b","-o", read_seq])

					print("\t".join([last_id, gene_folder, gene_folder+last_id,
									 left_output_seq, right_output_seq, read_seq]),
						  file=result_file_fd)

				gene_folder=result_folder + gene_id + '/'
				if ( not os.path.exists(gene_folder)):
					os.makedirs(gene_folder)
				tranfile=open(gene_folder+gene_id,"w+")
				last_id=gene_id
				chr=info[0]
				left_region=int(info[3])
				right_region=int(info[4])


		info[6]='+'
		print("\t".join(info)+"\n",file=tranfile)
	return 


def main():
	try:
		
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
		
		args=parser.parse_args()
		
		if (args.json != None and os.path.exists(args.json)):
			config_data   = open('config.json')
			configuration = json.load(config_data)
			result_folder	= configuration['result_folder'] + '/' + configuration['alias'] + '/'
			output_info_file = result_folder + configuration['alias']
			left_strain_ref_seq = configuration['left_strain_ref_seq']
			right_strain_ref_seq = configuration['right_strain_ref_seq']
			left_strain_id = configuration['left_strain_id']
			right_strain_id = configuration['right_strain_id']
			bam_file = configuration['bam_file']
			result_file	 =  result_folder+configuration['result_file']
			split_gtf_file(left_strain_ref_seq, right_strain_ref_seq,
						   left_strain_id, right_strain_id,
						   bam_file,
						   result_folder, output_info_file, result_file)
			return

		if (args.filelist != None and os.path.exists(args.filelist)):
			# if (args.reftable == None):
			# 	print('Error: no reftable supplied (--reftable)')
			# 	return 
				
			files=open(args.filelist).readlines()
			resultfolder='result/'+os.path.basename(args.filelist)+'/cufflinks/'
			if not os.path.exists(resultfolder):
				os.makedirs(resultfolder)
			for f in files:
				resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
				if not os.path.exists(resultsubfolder):
					os.makedirs(resultsubfolder)
				print('cufflinks -o '+resultsubfolder+' '+f.strip())
		
	except: 
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		print("*** print_exc:")
		traceback.print_exc()

if __name__ == "__main__":
	main()
		
