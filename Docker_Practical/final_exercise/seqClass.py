#!/usr/bin/env python

import sys, re
from argparse import ArgumentParser

parser = ArgumentParser(description = 'Classify a sequence as DNA or RNA')
parser.add_argument("-s", "--seq", type = str, required = True, help = "Input sequence")
parser.add_argument("-m", "--motif", type = str, required = False, help = "Motif")

# to show help if no argument was provided
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

args.seq = args.seq.upper()                 


if re.search('^[ACGTU]+$', args.seq):   #check that sequence contains only valid nucleotides

    if re.search('T', args.seq) and re.search('U', args.seq):   #check if the seq contains both T and U 
        print('The sequence is not valid DNA nor RNA (contains both T and U)')

    elif re.search('T', args.seq):   # if only T is present
        print ('The sequence is DNA')

    elif re.search('U', args.seq):   # if only U is oresent
        print ('The sequence is RNA')

    else:                            # if no T and no U present
        print ('The sequence can be DNA or RNA')
else:
    print ('The sequence is not DNA nor RNA')


# to check the motif
if args.motif:
    args.motif = args.motif.upper()
    print(f'Motif search enabled: looking for motif "{args.motif}" in sequence "{args.seq}"... ', end = '')
    if motif in args.seq:
        print("FOUND")
    else:
        print("NOT FOUND")
