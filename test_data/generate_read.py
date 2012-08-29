
import os
import json
import sys
import traceback
import argparse
import random
from subprocess import call


fa = open("test_transcript_seq.fa").readlines()
seq = fa[1]

left = open("left.fq", "w+")
right = open("right.fq", "w+")

for i in xrange(10000):
    length = random.randint(40,70)
    a = random.randint(0, len(seq) - length - 1)
    b = random.randint(a, len(seq) - length - 1)
    left.write(">read%s\n" % i)
    left.write(seq[a:(a + length)] + "\n")
    right.write(">read%s\n" % i)
    right.write(seq[b:(b + length)] + "\n")
