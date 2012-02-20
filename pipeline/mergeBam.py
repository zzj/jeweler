##copyright: kemal

import pysam
import sys
import numpy
import time
import os, os.path
import random

VERSION = "0.2"
PROPER_READS = (81, 83, 97, 99, 145, 147, 161, 163)
WEIGHT = (2,8,64)
CHUNKSIZE = 10000000
part = [dict() for i in xrange(3)]

def encodeQname(rseq):
    k = rseq.qname.rfind('#')
    name = rseq.qname[0:k]
    field = name.split(':')
    w = 0
    for i in xrange(3):
        if field[i] not in part[i]:
            k = len(part[i])
            part[i][field[i]] = k 
        else:
            k = part[i][field[i]]
        w += WEIGHT[i]*k
    w += 0 if rseq.is_read1 else 1
    x = int(field[3])
    y = int(field[4])
    # w range = 2 reads*4 machines*8 lanes*48 tiles = 0..3072  (12 bits)
    # x range = 0..21000    (15 bits)
    # y range = 0..200000   (18 bits)
    key = (w << 33) + (x << 18) + y
    # rseq.tid range = 0..100          (7 bits)
    # rseq.pos range = 0..200000000    (28 bits)
    value = (rseq.tid << 28) + rseq.pos
    return key, value

def usage(program):
    print ("Usage:")
    print ("    ",program," [-?] [-sort] [-v] infileA infileB [outfile]")
    print ("-sort: do not sort final output")
    print ("-v: output verbose messages")
    print ("-?: Output this usage guide")
    
def parseCommandLine(args):
    filenameA = ''
    filenameB = ''
    outfile = ''
    flags = {'verbose':False, 'sort':True}
    i = 1
    while (i < len(args)):
        arg = args[i]
        if (arg == '-?') or (arg == '--help'):
            usage(args[0])
            sys.exit(0)
        elif (arg == '-sort'):
            flags['sort'] = False
        elif (arg == '-v'):
            flags['verbose'] = True
        elif (filenameA == ''):
            filenameA = arg
        elif (filenameB == ''):
            filenameB = arg
        elif (outfile == ''):
            outfile = arg
        else:
            usage(args[0])
            sys.exit(1)
        i += 1
    if (filenameA == '') or (filenameB == ''):
        print ("ERROR: missing input files")
        usage(args[0])
        sys.exit(1)
    if (outfile == ''):
        fA = os.path.basename(filenameA)
        fB = os.path.basename(filenameB)
        outfile = os.path.commonprefix([fA,fB]) + 'Merged.bam'
    return filenameA, filenameB, outfile, flags

def mergeBams(args):
    filenameA, filenameB, outfile, flag = parseCommandLine(args)
    # Pass 1: scan through file B, make a cache of the reads it contains,
    #         and where they map to
    if flag['verbose']:
        start = time.time()
        print ("Starting Pass #1: filenameB")
    readCache = numpy.empty((2,CHUNKSIZE), dtype='int64')
    nB = 0
    bamB = pysam.Samfile(filenameB, 'rb')
    for rseq in bamB:
        key, value = encodeQname(rseq)
        readCache[:,nB] = value, key
        nB += 1
        if (nB == readCache.shape[1]):
            expand = numpy.empty((2,CHUNKSIZE), dtype='int64')
            readCache = numpy.hstack((readCache,expand))
    bamB.close()
    order = numpy.lexsort(readCache[:,0:nB])
    readCache = readCache[:,order]
    
    # Pass 2: scan through file A, output every reads it contains while
    #         checking if the same read is also mapped in file B
    if flag['verbose']:
        t = time.time()
        print( "Pass 1: ",nB," reads (",t-start," secs)" )
        start = t
        print( "Starting Pass #2: ", filenameA)
    nA = 0
    nU = 0
    setFound = 1 << 36
    bamA = pysam.Samfile(filenameA, 'rb')
    outHeader = dict(bamA.header.items())
    outHeader['PG'] = [{'ID': 'annotateSNPs', 'PP': outHeader['PG'][0]['ID'], 'VN': VERSION, 'CL': ' '.join(args)}] + outHeader['PG']
    tmpFilename = "/tmp/mrg%08d.bam" % random.randint(0, 99999999)
    bamO = pysam.Samfile(tmpFilename, 'wb', header=outHeader, referencenames=bamA.references)
    for rseq in bamA:
        key, value = encodeQname(rseq)
        i = numpy.searchsorted(readCache[1,:], (key), side='left')
        if (i < nB) and (readCache[1,i] == key):
            j = numpy.searchsorted(readCache[1,:], (key), side='right')
            for k in xrange(i,j):
                if (readCache[0,k] == value):
                    # Read maps to the same location in both alignments A and B
                    rseq.tags += [('YA',3)]
                    readCache[0,k] += setFound
                    nU += 1
                    break
            else:
                # Note: this is a case needs to be handled better.
                # This read appears in file B, but maps to a unique
                # position in file A. If it mapped to a single postion
                # in file B (NH = 1), then we are creating a repeat.
                # If it is repeated in file B (NH > 1), then we should
                # increase the NH flag on this and the all other copies of
                # this read.
                #
                # The fix might require waiting to the 3rd pass to handle
                # this case.
                rseq.tags += [('YA',1)]
        else:
            # Read is unique to alignment A
            rseq.tags += [('YA',1)]
        bamO.write(rseq)
        nA += 1
    bamA.close()

    # Pass 3: scan through file B, output those reads that are unique to it
    #
    if flag['verbose']:
        t = time.time()
        print ("Pass 2: ",nA," reads, %d in intersection (",t-start," secs)")
        start = t
        print ("Starting Pass #3: ",  outfile)
    nB = 0
    bamB = pysam.Samfile(filenameB, 'rb')
    for rseq in bamB:
        key, value = encodeQname(rseq)
        i = numpy.searchsorted(readCache[1,:], (key), side='left')
        j = numpy.searchsorted(readCache[1,:], (key), side='right')
        for k in xrange(i,j):
            if (readCache[0,k] == value):
                # Read is unique to alignment B, if it was also
                # in A than setFound would be added to value and
                # it would not have matched
                rseq.tags += [('YA',2)]
                nB += 1
                bamO.write(rseq)
                break
    bamB.close()
    bamO.close()

    # A final sorting step
    #
    if flag['verbose']:
        t = time.time()
        print( "Pass 3: ",nB," reads (",t-start," secs)" )
    if flag['sort']:
        if flag['verbose']:
            start = t
            print ("Sorting result")
        pysam.sort(tmpFilename, outfile[:outfile.rfind('.bam')])
        os.remove(tmpFilename)
        if flag['verbose']:
            t = time.time()
            print ("Sort: (",t-start," secs)")
    else:
        os.rename(tmpFilename, outfile)

if __name__ == '__main__':
    mergeBams(sys.argv)