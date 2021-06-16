##Libraries
import argparse
import numpy as np
import copy
from datetime import datetime
from collections import defaultdict
import csv
import pandas as pd

#Define arguments for each required and optional input
def defineArguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--gb-processed-file",  dest="GbProcessedFile",required=True,help="GbProcessedFile")
    parser.add_argument("--clustering-only",  dest="ClusterOnly",required=False,action="store_true",help="ClusterOnly")
    parser.add_argument("--kmer-size", dest="KmerSize",required=False,default=10,help="KmerSize",type=int)
    parser.add_argument("--num-kmers-per-sequence", dest="NumKmersPerSequence",required=False,default=10,help="NumKmersPerSequence",type=int)
    parser.add_argument("--num-clusters", dest="NumClusters",required=False,default=10,help="NumClusters",type=int)
    parser.add_argument("--output-directory",  dest="OutputDirectory",required=True,help="OutputDirectory")

    return parser

def parse_files(GbProcessedFile):

    gbData = pd.read_csv(GbProcessedFile, sep = '\t')

    return gbData

def declare_AA_alphabet():

    AAhash = { 'T': 'A',
          'S': 'A',
          'R': 'K',
          'D': 'N',
          'Q': 'E',
          'V': 'I',
          'M': 'L',
          'Y': 'F'}

    return AAhash

def replace_aas(AAhash,gbData):

    for i in range(len(gbData)):
        for j in (AAhash.keys()):
            gbData['AAsequence'][i] = gbData['AAsequence'][i].replace(j,AAhash[j])

    return gbData

def main():

    #Generate argument parser and define arguments
    parser = defineArguments()
    args = parser.parse_args()
    
    GbProcessedFile = args.GbProcessedFile
    ClusterOnly = args.ClusterOnly
    KmerSize = args.KmerSize
    NumKmersPerSequence = args.NumKmersPerSequence
    NumClusters = args.NumClusters
    OutputDirectory = args.OutputDirectory

    #Command expected to execute the script
    cmd = "--gb-processed-file %s --kmer-size %i --num-kmers-per-sequence %i --num-clusters %i --output-directory %s" % (GbProcessedFile,KmerSize,NumKmersPerSequence,NumClusters,OutputDirectory)
    



if __name__ == '__main__':
    main()
