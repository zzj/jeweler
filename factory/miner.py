import os
import json
import sys
import traceback
import argparse
from subprocess import call

from  split_gtf_file import split_gtf_file
from  miner_parser import initialize_parser

def main():
	try:
		parser = initialize_parser()

		args=parser.parse_args();

		if (args.is_single):
			## process a single bam file
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
					## run cufflinks command
					resultfolder='result/'+os.path.basename(args.filelist)+'/cufflinks/'
					if not os.path.exists(resultfolder):
						os.makedirs(resultfolder)
					for f in files:
						resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
						if not os.path.exists(resultsubfolder):
							os.makedirs(resultsubfolder)
						print('cufflinks -o '+resultsubfolder+' '+f.strip())
					
				elif (args.is_jeweler):
					resultfolder = 'result/'+os.path.basename(args.filelist)+'/jeweler/'
					appraiser_resultfolder = 'result/'+os.path.basename(args.filelist)+'/appraiser/'
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
						if (args.is_jeweler_only):
							## multiple alignments maps
							mamf_command = ""
							if (os.path.exists(appraiser_resultfolder + alias +"/sm_table")):
								mamf_command = " -mamf " + appraiser_resultfolder + alias +"/sm_table";

							print('./jeweler -i '+result_folder+'/'+alias+'/'+result_file + mamf_command + " -earrings")
						elif (args.is_mismatch_analyzer):
							ma_resultfolder = 'result/'+os.path.basename(args.filelist)+'/mismatch_analyzer/'
							if not os.path.exists(ma_resultfolder):
								os.makedirs(ma_resultfolder)
							ma_resultfolder = 'result/'+os.path.basename(args.filelist)+'/mismatch_analyzer/' + alias + '/'
							if not os.path.exists(ma_resultfolder):
								os.makedirs(ma_resultfolder)
							print('./jeweler -i '+result_folder+'/'+alias+'/'+result_file +  " -mismatch_analyzer " + ma_resultfolder+ "result")
						else:
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
				elif (args.is_appraiser) :
					resultfolder='result/'+os.path.basename(args.filelist)+'/appraiser/'
					if not os.path.exists(resultfolder):
						os.makedirs(resultfolder)
					for f in files:
						resultsubfolder=resultfolder+'/'+os.path.basename(f.strip().replace('.bam',''))
						if not os.path.exists(resultsubfolder):
							os.makedirs(resultsubfolder)
						print('./appraiser -mamf '+resultsubfolder+'/sm -l '+resultsubfolder+'/log -qf '+ resultsubfolder +'/qf -b '+f.strip())
					

	except: 
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		print("*** print_exc:")
		traceback.print_exc()

if __name__ == "__main__":
	main()
		
