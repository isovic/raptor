#! /usr/bin/env python2.7

import re
import sys
import argparse
import fastqparser

################################
######### Utility tools ########
import sys;
import os;
from time import gmtime, strftime
import traceback;
import fnmatch;

### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log=sys.stderr):
  timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());      # pragma: no cover
  if (fp_log != None):                                      # pragma: no cover
      fp_log.write('[%s] %s\n' % (timestamp, message))      # pragma: no cover
      fp_log.flush();                                       # pragma: no cover
################################

COMPLEMENT = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}

def revcmp(seq):
    return ''.join([COMPLEMENT[v] for v in seq[::-1]])

def fastq_invert(reads, region_start, region_end, fp_out):
    s = region_start if region_start >= 0 else 0

    for seq_fields in fastqparser.yield_seq([reads]):
        header, seq = seq_fields[0:2]
        e = region_end if region_end < len(seq) else len(seq)

        prefix = ''
        infix = ''
        suffix = ''

        if s < len(seq):
            prefix = seq[:s]
            infix = revcmp(seq[s:e])
            new_header = '%s-inverted_%d-%d' % (header, s, e)
        else:
            prefix = seq
            infix = ''
            new_header = '%s-not_inverted' % (header)
        if e < len(seq):
            suffix = seq[e:]

        new_seq = prefix + infix + suffix

        fp_out.write('%s\n%s\n' % (new_header, new_seq))

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Inverts a specified region to all reads in the specified input dataset, if the read is long enough to have this region.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('region', help='Region specified in the format: "start-end".')
    parser.add_argument('reads', help='A set of reads in FASTA/FASTQ format.')
    parser.add_argument('--out', type=str, default='-', help='Modified reads.')

    args = parser.parse_args(argv[1:])

    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    assert(args.out != args.reads)    # Prevent accidental overwriting.

    fp_out = sys.stdout if (args.out == '-') else open(args.out, 'w')

    region_start, region_end = args.region.split('-')
    region_start = int(region_start)
    region_end = int(region_end)

    fastq_invert(args.reads, region_start, region_end, fp_out)

    if fp_out != None and fp_out != sys.stdout: # pragma: no cover
        fp_out.close()                          # pragma: no cover

if __name__ == "__main__":  # pragma: no cover
    main(sys.argv)          # pragma: no cover
