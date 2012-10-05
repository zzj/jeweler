import os
def split_fasta(filename, folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
    fd = None
    for line in open(filename, 'r'):
        if line[0] == ">":
            if fd != None:
                fd.close()
                fd = None
            fd = open(folder + "/" + line.strip()[1:] + ".fa", "w")
        fd.write(line)
if __name__ == "__main__":
    split_fasta("../data/index/CAST.fa", "../data/index/CAST/")
    split_fasta("../data/index/WSB.fa", "../data/index/WSB/")
    split_fasta("../data/index/PWK.fa", "../data/index/PWK/")