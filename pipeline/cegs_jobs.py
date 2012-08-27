import os, sys, socket, time

def run_command(fd, comm):
    print(comm)
    ##os.system(comm)

def main():
    job_id = int(sys.argv[1])
    log_file_fd = open('log/log_'+str(job_id),"w+")
    tophat_log_file = 'log_tophat'+str(job_id)

    playpen_index = '/lustre/scr/z/z/zzj/jeweler/data/index/'
    playpen_root = '/nas02/home/z/z/zzj/Research/rna_seq/jeweler/data/cegs_rnaseq_bam_new/'
    fastq_root = '/lustre/scr/k/e/kemal/RNAseq/data/'
    ##run_command(log_file_fd,'mkdir '+playpen_root)

    ref_map={'F':'CAST', 'G':'PWK', 'H':'WSB'}

    with open('../info/fastq') as f:
        read_meta = f.readlines()

    first_read_file = read_meta[job_id * 2 - 2].strip()
    second_read_file = read_meta[job_id * 2 - 1].strip()
    first_local_read_file = fastq_root + first_read_file
    second_local_read_file = fastq_root + second_read_file

    filebasename = os.path.basename(first_read_file)
    batch = first_read_file.strip().split('_')[2].split('/')[0]

    info = filebasename.strip().split('_')
    lane_id = info[info.index('s') + 1]

    ##print(lane_id + '\n')

    output_folder = info[0][0:2] + '_' + info[0][2:8] + '_' + info[1] + '_aligned_to_' + ref_map[info[0][0]] + '_lane' + lane_id + '_' + batch
    local_output_folder = playpen_root + output_folder
    run_command(log_file_fd,'tophat  -o ' + local_output_folder + ' -N 8 -p 4 -r 100 ' + playpen_index+ref_map[info[0][0]] + '.fa ' + first_local_read_file + ' ' + second_local_read_file)

    if (info[0][0]!=info[0][1]):
        output_folder=info[0][0:2]+'_'+info[0][2:8]+'_'+info[1]+'_aligned_to_'+ref_map[info[0][1]]+'_lane'+lane_id+'_'+batch
        local_output_folder = playpen_root + output_folder
        run_command(log_file_fd, 'tophat  -o ' + local_output_folder + ' -N 8 -p 4 -r 100 ' + playpen_index + ref_map[info[0][1]] + '.fa ' + first_local_read_file + ' ' + second_local_read_file)

    log_file_fd.close()


if __name__ == "__main__":
    main()
