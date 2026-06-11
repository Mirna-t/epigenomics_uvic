#!/usr/bin/env python

import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)

x = 1000000
selectedGene = ""
selectedGeneStart = 0

for line in open_input.readlines():
	gene, y = line.strip().split('\t')

	position = int(y)
	distance = abs(position - enhancer_start)

	if distance < x:
		x = distance
		selectedGene = gene
		selectedGeneStart = position

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])

