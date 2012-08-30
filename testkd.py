import os
import sys

QUEUE = "day"
MEM = "8"
PROCESS = "1"
ID = 1
def pause():
    while True:
        p = os.popen("bjobs")
        output = p.readlines()
        p.close()
        if len(output) <= 1:
            return
        else:
            print("sleeping")
            os.system("sleep 100")
            continue


def config(f):
    global QUEUE, MEM, PROCESS
    _, QUEUE, MEM, PROCESS = f.split(" ")


def run_command(f):
    global ID
    cmd = "bsub -C 0 -q " + QUEUE + " -M " + MEM + " -n " + PROCESS + " -R \"span[hosts=1]\"  -o temp/test" + str(ID) + " " + f
    # print (cmd)
    os.system(cmd)
    ID += 1


def process_line(f):
    if f == "pause":
        pause()
    elif f[0:6] == "config":
        config(f)
    else:
        run_command(f)


def main():
    files = open("command_list").readlines()
    for f in files:
        process_line(f.strip())
        
if __name__ == "__main__":
    main()
