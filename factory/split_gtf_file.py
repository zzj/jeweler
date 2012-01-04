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
					maternal_output_seq  = gene_folder+last_id+'.'+maternal_strain_id+".fa"
					paternal_output_seq = gene_folder+last_id+'.'+paternal_strain_id+".fa"
					read_seq		 = gene_folder+last_id+".seq.bam";
					call(['gffread', '-w', maternal_output_seq, '-g', maternal_strain_ref_seq, 
						  gene_folder+last_id])
					call(['gffread', '-w', paternal_output_seq, '-g', paternal_strain_ref_seq, 
						  gene_folder+last_id])
					call(["samtools","view", bam_file , 
						 chr +":" +str(left_region)  +"-" +str(right_region) ,
						 "-b","-o", read_seq])

					print("\t".join([last_id, gene_folder, gene_folder+last_id,
									 maternal_output_seq, paternal_output_seq, read_seq]),
						  file=result_file_fd)
					result_file_fd.flush();

				gene_folder=result_folder + gene_id + '/'
				if ( not os.path.exists(gene_folder)):
					os.makedirs(gene_folder)
				tranfile=open(gene_folder+gene_id,"w+")
				last_id=gene_id
				chr=info[0]
				left_region=int(info[3])
				right_region=int(info[4])


		info[6]='+'
		print("\t".join(info),file=tranfile)
	result_file_fd.close()
	##os.system('./jeweler -i '+result_file)
	return 

