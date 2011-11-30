import os
import json
import sys
import traceback
from subprocess import call
try:
    config_data          = open('config.json')
    configuration        = json.load(config_data)
    result_folder        = configuration['result_folder'] + '/' + configuration['alias'] + '/'
    output_info_file     = result_folder + configuration['alias']
    left_strain_ref_seq  = configuration['left_strain_ref_seq']
    right_strain_ref_seq = configuration['right_strain_ref_seq']
    left_strain_id       = configuration['left_strain_id']
    right_strain_id      = configuration['right_strain_id']
    bam_file             = configuration['bam_file']
    result_file          = result_folder+configuration['result_file']  
    result_file_fd       = open(result_file,"w+")
    try:
        os.makedirs()
    except:
        print("Notice: the folder was created.\n")


    with open(configuration['gtf_input_file']) as f:
        gene_meta = f.readlines()

    last_id          = ""
    chr             = ""
    left_region     = -1
    right_region    = -1
    is_file_created = False
    tranfile        = None
    id              = 0

    for line in gene_meta:
        info=line.strip(' \t\n\r').split('\t')
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
                    read_seq         = gene_folder+last_id+".seq.bam";
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

except: 
    exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
    print("*** print_exc:")
    traceback.print_exc()
