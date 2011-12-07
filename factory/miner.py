import os
import json
import sys
import traceback
import argparse
from subprocess import call


def split_gtf_file(gtf_input_file,maternal_strain_ref_seq, paternal_strain_ref_seq,
					maternal_strain_id, paternal_strain_id,
					bam_file,
					result_folder,  result_file):
	result_file_fd= open(result_file,"w+")
	try:
		os.makedirs()
	except:
		print("Notice: the folder was created.\n")


	with open(gtf_input_file) as f:
		gene_meta = f.readlines()

	last_id		  = ""
	chr			 = ""
	maternal_region	 = -1
	paternal_region	= -1
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
				if ((maternal_region)<0):
					maternal_region = int(info[3])
				else :
					maternal_region = min(int(info[3]),maternal_region)

				if (paternal_region<0):
					paternal_region=int(info[4])
				else :
					paternal_region=max(int(info[4]),paternal_region)
			else :
				# new transcript, new gene
				if (tranfile!=None):
					tranfile.close()
					id=id+1
					print(str(id)+"\n")
					maternal_output_seq  = gene_folder+last_id+'.'+maternal_strain_id+".fa"
					paternal_output_seq = gene_folder+last_id+'.'+paternal_strain_id+".fa"
					read_seq		 = gene_folder+last_id+".seq.bam";
					call(['gffread', '-w', maternal_output_seq, '-g', maternal_strain_ref_seq, 
						  gene_folder+last_id])
					call(['gffread', '-w', paternal_output_seq, '-g', paternal_strain_ref_seq, 
						  gene_folder+last_id])
					call(["samtools","view", bam_file , 
						 chr +":" +str(maternal_region)  +"-" +str(paternal_region) ,
						 "-b","-o", read_seq])

					print("\t".join([last_id, gene_folder, gene_folder+last_id,
									 maternal_output_seq, paternal_output_seq, read_seq]),
						  file=result_file_fd)

				gene_folder=result_folder + gene_id + '/'
				if ( not os.path.exists(gene_folder)):
					os.makedirs(gene_folder)
				tranfile=open(gene_folder+gene_id,"w+")
				last_id=gene_id
				chr=info[0]
				maternal_region=int(info[3])
				paternal_region=int(info[4])


		info[6]='+'
		print("\t".join(info),file=tranfile)
	result_file_fd.close()
	os.system('./jeweler -i '+result_file)
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

		parser.add_argument('--single',
							action='store_true',
							dest='is_single')

		parser.add_argument('--cufflinks',
							action='store_true',
							dest='is_cufflinks')

		parser.add_argument('--jeweler',
							action='store_true',
							dest='is_jeweler')



		args=parser.parse_args();

		if (args.is_single):
			configuration=dict()
			if (args.json != None and os.path.exists(args.json)):
				config_data   = open(args.json)
				configuration = json.load(config_data)
			else :
				if (args.maternal_strain_id is None or
					args.paternal_strain_id is None or
					args.maternal_strain_ref_seq is None or
					args.paternal_strain_ref_seq is None or
					args.result_file is None or
					args.result_folder is None or
					args.alias is None or
					args.gtf_input_file is None):
					print("Missing arguments")
					return
				configuration['maternal_strain_ref_seq']=args.maternal_strain_ref_seq
				configuration['maternal_strain_id']=args.maternal_strain_id
				configuration['paternal_strain_ref_seq']=args.paternal_strain_ref_seq
				configuration['paternal_strain_id']=args.paternal_strain_id
				configuration['result_file']=args.result_file
				configuration['alias']=args.alias
				configuration['bam_file']=args.bam_file
				configuration['result_folder']=args.result_folder
				configuration['gtf_input_file']=args.gtf_input_file
				
			result_folder	= configuration['result_folder'] + '/' + configuration['alias'] + '/'
			maternal_strain_ref_seq = configuration['maternal_strain_ref_seq']
			paternal_strain_ref_seq = configuration['paternal_strain_ref_seq']
			maternal_strain_id = configuration['maternal_strain_id']
			paternal_strain_id = configuration['paternal_strain_id']
			bam_file = configuration['bam_file']
			result_file	 =  result_folder+configuration['result_file']
			gtf_input_file=configuration['gtf_input_file']
			split_gtf_file(gtf_input_file, maternal_strain_ref_seq, paternal_strain_ref_seq,
						   maternal_strain_id, paternal_strain_id,
						   bam_file,
						   result_folder, result_file)
			return
				
		else :

			refidtable=dict()
			reffiletable=dict()
			if (args.reftable is not None and os.path.exists(args.reftable)):
				with open(args.reftable) as f:
					refs = f.readlines()
				for ref in refs:
					r=ref.strip().split('\t')
					refidtable[r[0]]=r[1]
					reffiletable[r[0]]=r[2]

			else :
				print('No reftable supplied')
				return 
			if (args.filelist != None and os.path.exists(args.filelist)):
				# if (args.reftable == None):
				# 	print('Error: no reftable supplied (--reftable)')
				# 	return 

				files=open(args.filelist).readlines()
				if (args.is_cufflinks):
					resultfolder='result/'+os.path.basename(args.filelist)+'/cufflinks/'
					if not os.path.exists(resultfolder):
						os.makedirs(resultfolder)
					for f in files:
						resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
						if not os.path.exists(resultsubfolder):
							os.makedirs(resultsubfolder)
						print('cufflinks -o '+resultsubfolder+' '+f.strip())
				elif (args.is_jeweler):
					resultfolder='result/'+os.path.basename(args.filelist)+'/jeweler/'
					if not os.path.exists(resultfolder):
						os.makedirs(resultfolder)
					for f in files:
						resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
						if not os.path.exists(resultsubfolder):
							os.makedirs(resultsubfolder)
						result_folder=resultfolder
						alias=os.path.basename(f.strip().replace('.bam',''))
						maternal_strain_id=refidtable[alias[0]]
						paternal_strain_id=refidtable[alias[1]]
						maternal_strain_ref_seq=reffiletable[alias[0]]
						paternal_strain_ref_seq=reffiletable[alias[1]]
						bam_file=f.strip()
						result_file=alias+'.info'
						cuffresultfolder='result/'+os.path.basename(args.filelist)+'/cufflinks/'
						gtf_input_file=	cuffresultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))+'/transcripts.gtf'
						print('python3.2 factory/miner.py --single '+
							  '--maternal_strain_ref_seq '+maternal_strain_ref_seq +' ' +
							  '--paternal_strain_ref_seq '+paternal_strain_ref_seq +' ' +
							  '--maternal_strain_id '+maternal_strain_id +' ' +
							  '--paternal_strain_id '+paternal_strain_id +' ' +
							  '--bam_file '+bam_file +' ' +
							  '--alias '+alias +' ' +
							  '--result_file '+result_file +' ' +
							  '--result_folder '+result_folder +' '+
							  '--gtf_input_file '+gtf_input_file
							)

					
					
		
		
	except: 
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		print("*** print_exc:")
		traceback.print_exc()

if __name__ == "__main__":
	main()
		
