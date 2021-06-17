##Libraries
import argparse
import numpy as np
import copy
from datetime import datetime
from collections import defaultdict
import csv
import pandas as pd
import heapq
import matplotlib.pyplot as plt

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

def compute_kmer_hashes(gbData,KmerSize,NumKmersPerSequence):

    top_kmer_dict = []

    #iterate through all values in gbData
    for i in range(len(gbData)):
    #for i in range(1):
        all_kmer_dict = []
        #Calc hash for the first kmer found in the sequence
        new_hash = {}
        kmer_hash = 0
        AAsequence = gbData['AAsequence'][i]
        kmer = AAsequence[0:KmerSize]
        for j in range(0,KmerSize-1):
            kmer_hash += ord(kmer[j])**(KmerSize - (j+1))
        kmer_hash += ord(kmer[KmerSize-1])
        new_hash['Kmer_hash'] = kmer_hash
        new_hash['Kmer'] = kmer
        new_hash['Locus'] = gbData['LOCUS'][i]
        new_hash['Index'] = 0
        new_hash['Seq_length'] = len(AAsequence)
        all_kmer_dict.append(new_hash)
        new_hash = {}
        
        #Calc hash for all other kmers found in the sequence
        for m in range(1,len(AAsequence) - KmerSize + 1):
            kmer = AAsequence[m:m+KmerSize]
            kmer_hash = kmer_hash - (ord(AAsequence[m-1])**(KmerSize-1)) + ord(kmer[KmerSize-1])
            new_hash['Kmer_hash'] = kmer_hash
            new_hash['Kmer'] = kmer
            new_hash['Locus'] = gbData['LOCUS'][i]
            new_hash['Index'] = m
            new_hash['Seq_length'] = len(AAsequence)
            all_kmer_dict.append(new_hash)
            new_hash = {}
            
        #Subset Kmer Hash based on top M hash vals
        hash_vals = [x['Kmer_hash'] for x in all_kmer_dict]
        top_vals = heapq.nlargest(NumKmersPerSequence,hash_vals)
        all_kmer_dict = [d for d in all_kmer_dict if d['Kmer_hash'] in top_vals]
            
        top_kmer_dict = top_kmer_dict + all_kmer_dict

    return top_kmer_dict

def assign_clusters(top_kmer_dict):

    unique_hash = []

    for i in range(len(top_kmer_dict)):
        if(top_kmer_dict[i]['Kmer_hash'] not in unique_hash):
            unique_hash.append(top_kmer_dict[i]['Kmer_hash'])
            top_kmer_dict[i]['Cluster'] = (unique_hash.index(top_kmer_dict[i]['Kmer_hash'])+1)
        else:
            top_kmer_dict[i]['Cluster'] = (unique_hash.index(top_kmer_dict[i]['Kmer_hash'])+1)

    return unique_hash

def find_centers(top_kmer_dict):

    all_clusters = [x['Cluster'] for x in top_kmer_dict]
    unique_clusters = set(all_clusters)

    cluster_centers = {}

    for i in unique_clusters:
        #val = {}
        cluster = [i]
        clust_kmer_dict = [d for d in top_kmer_dict if d['Cluster'] in cluster]
        seq_lengths = [x['Seq_length'] for x in clust_kmer_dict] 
        top_length = max(seq_lengths)
        
        #Get first seq length that matches the top one
        for j in range(len(clust_kmer_dict)):
            if(clust_kmer_dict[j]['Seq_length'] == top_length):
                cluster_centers[clust_kmer_dict[j]['Cluster']] = clust_kmer_dict[j]['Locus']
                #Only use first sequence of that length
                break

    #Assign new cluster to each unique center
    all_centers = pd.DataFrame(cluster_centers.values())
    all_centers = all_centers.drop_duplicates()
    all_centers = all_centers.reset_index(drop=True)
    combined_clusters = {}

    clust = 1
    for i in all_centers[0]:
        combined_clusters[i] = clust
        clust += 1

    return cluster_centers,combined_clusters

def reassign_clusters(top_kmer_dict,cluster_centers,combined_clusters):

    for i in range(len(top_kmer_dict)):
        
        #If locus is not center locus for this cluster
        if(top_kmer_dict[i]['Locus'] != cluster_centers[top_kmer_dict[i]['Cluster']]):
            #Which locus does the cluster belong to
            center_locus = cluster_centers[top_kmer_dict[i]['Cluster']]
            top_kmer_dict[i]['Cluster'] = combined_clusters[center_locus]
        else:
            top_kmer_dict[i]['Cluster'] = combined_clusters[top_kmer_dict[i]['Locus']]
    
    #Convert to dataframe with only locus and cluster assignment
    protein_clusters_df = pd.DataFrame(pd.DataFrame(top_kmer_dict), columns = ['Locus','Cluster'])
    protein_clusters_df = protein_clusters_df.drop_duplicates(subset=['Locus'])
    protein_clusters_df = protein_clusters_df.reset_index(drop=True)

    #Add metadata
    protein_clusters_meta_df = pd.concat([protein_clusters_df,gbData[['Location', 'collectionDate']].reset_index(drop=True)],axis=1)

    protein_clusters_df.to_csv(OutputDirectory+"/Clusters.csv", index = False)
    protein_clusters_meta_df.to_csv(OutputDirectory+"/Clusters_Meta.csv", index = False)

    return protein_clusters_df,protein_clusters_meta_df

def output_cluster_file(cluster_file,OutputDirectory):

    cluster_file.to_csv(OutputDirectory+"/Clusters.csv", index = False)

def cluster_histogram(protein_clusters_df,OutputDirectory):

    x = list(protein_clusters_df['Cluster'])

    plt.hist(list(protein_clusters_df['Cluster']),bins = len(set(x)))
    axes = plt.gca()
    axes.set_ylim([0,1500])
    #Save figure to user specified output folder.
    plt.savefig(OutputDirectory+"/ClusterHistogram.png",bbox_inches='tight')
    plt.clf()

##Create files for Cytoscape
def create_cluster_edge_list(combined_clusters,protein_clusters_df):

    clusters_edge_list = pd.DataFrame()

    for i in list(combined_clusters.values()):
        clusters_set = protein_clusters_df.loc[protein_clusters_df['Cluster'] == i]
        for j in range(len(clusters_set)):
            list_range = [item for item in range(j+1, len(clusters_set))]
            for k in list_range:
                clusters_edge_list = clusters_edge_list.append({'Source':clusters_set['Locus'].iloc[j],'Target':clusters_set['Locus'].iloc[k]},ignore_index=True)
    
    return clusters_edge_list


def output_cytoscape_files(clusters_edge_list,cluster_file,OutputDirectory):

    clusters_edge_list.to_csv(OutputDirectory+"/Clusters_EdgeList.csv", index = False)
    cluster_file.to_csv(OutputDirectory+"/Clusters_Metadata.noa", index = False)



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

    starttime = datetime.now()
    print("Sars_Protein_Clustering Start Time: ",starttime)

    gbData = parse_files(GbProcessedFile)

    #Prepare sequences
    AAhash = declare_AA_alphabet()
    gbData = replace_aas(AAhash,gbData)

    #Compute hash values for each kmer
    top_kmer_dict = compute_kmer_hashes(gbData,KmerSize,NumKmersPerSequence)

    print("Hashes computed: ",datetime.now())

    #Assign Clusters
    unique_hash = assign_clusters(top_kmer_dict)
    cluster_centers,combined_clusters = find_centers(top_kmer_dict)
    protein_clusters_df,protein_clusters_meta_df = reassign_clusters(top_kmer_dict,cluster_centers,combined_clusters)

    print("Clusters assigned: ",datetime.now())

    #Outputs
    output_cluster_file(protein_clusters_df,OutputDirectory)
    output_cluster_file(protein_clusters_meta_df,OutputDirectory)
    cluster_histogram(protein_clusters_df,OutputDirectory)

    #Cytoscape
    cluster_edge_list = create_cluster_edge_list(combined_clusters,protein_clusters_df)
    output_cytoscape_files(clusters_edge_list,cluster_file,OutputDirectory)

    endtime = datetime.now()
    tdelta = endtime-starttime
    print("Total computation time: ",tdelta)

if __name__ == '__main__':
    main()
