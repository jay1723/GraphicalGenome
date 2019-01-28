#! /usr/bin/python

import argparse
import CCGG

parser = argparse.ArgumentParser()

parser.add_argument("nodefile", help="Absolute path to the nodefile", type=str)
parser.add_argument("edgefile", help="Absolute path to the edgefile", type=str)
parser.add_argument("chromosome", help="The chromosome represented by the node and edge files", type=str)
parser.add_argument("--sequence", help="Extract the strain for the sequence. 2 required arguments: Strain, filepath", type=str, nargs=2)
args = parser.parse_args()


if args.extract:
    graph = CCGG.GraphicalGenome(args.nodefile, args.edgefile)
    strain = args.extract[0]
    filepath = args.extract[1]
    chromo = args.chromosome
    graph.writeSequence(strain, chromo, filepath) 