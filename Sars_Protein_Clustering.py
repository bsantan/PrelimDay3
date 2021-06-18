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
import seaborn as sns
import re

#Define arguments for each required and optional input
def defineArguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ##Need to add logic for how to handle these entries
    #When provided, will use output as gbData
    parser.add_argument("--gb-original-file",  dest="GbOriginalFile",required=False,help="GbOriginalFile")
    parser.add_argument("--output-processed-file",  dest="OutputProcessedFile",required=False,action="store_true",help="OutputProcessedFile")
    #When provided, will not perform parsing on original genBank file
    parser.add_argument("--gb-processed-file",  dest="GbProcessedFile",required=False,help="GbProcessedFile")
    parser.add_argument("--kmer-size", dest="KmerSize",required=False,default=10,help="KmerSize",type=int)
    parser.add_argument("--num-kmers-per-sequence", dest="NumKmersPerSequence",required=False,default=10,help="NumKmersPerSequence",type=int)
    parser.add_argument("--num-clusters", dest="NumClusters",required=False,default=10,help="NumClusters",type=int)
    parser.add_argument("--output-directory",  dest="OutputDirectory",required=True,help="OutputDirectory")

    return parser

def parse_original_file(GbOriginalFile):

    f = open(GbOriginalFile, "r")

    S_entries = []
    #Initialize the dict once
    entry_dict = {}
    desired_gene = False
    getting_protein = False
    seq = ''

    for i,line in enumerate(f):
        line = line.strip()
        values = line.split()
        if values:
            #Get all necessary items for this entry
            if values[0].lower() == 'accession':
                entry_dict['LOCUS'] = values[1]
            if 'country' in values[0].lower():
                #For USA locations with state listed
                if ":" in line:
                    values[0] = values[0]+values[1]
                #Get string inside quotation marks
                entry_dict['Location'] = re.findall('"([^"]*)"',values[0])[0]
            if 'collection_date' in values[0].lower():
                #Get string inside quotation marks
                entry_dict['collectionDate'] = re.findall('"([^"]*)"',values[0])[0]
            
            #Check if entry is of S protein and continue searching for next protein sequence
            if 'gene="s"' in values[0].lower():
                #Ensures only 1 protein sequence that matches S gene is found
                desired_gene = True
                
            if desired_gene is True:
                if 'translation' in values[0].lower():
                    #First check if protein is short and both double quotes are preset
                    if values[0].count('"') == 2:
                        #Get string inside quotation marks
                        entry_dict['AAsequence'] = re.findall('"([^"]*)"',values[0])[0]
                        #Only add entry if S protein was found
                        S_entries.append(entry_dict)
                        desired_gene = False
                        #Ensures that only lines that are part of protein sequence are appended
                        getting_protein = False
                        #Reset dict only after S protein has been found
                        entry_dict = {}
                        seq = ''
                    #Otherwise get first line of protein sequence
                    else:
                        seq = seq+line.split('"')[1]
                        print(seq)
                        getting_protein = True
                #Append all subsequence lines to protein sequence without ending quote
                if '"' not in line and getting_protein is True:
                    seq = seq+line
                #Append last line to protein sequence with ending quote
                if '"' in line and 'translation' not in line and getting_protein is True:
                    line = line.split('"')[0]
                    seq = seq+line
                    entry_dict['AAsequence'] = seq
                    #Only add entry if S protein was found
                    S_entries.append(entry_dict)
                    desired_gene = False
                    getting_protein = False
                    #Reset dict only after S protein has been found
                    entry_dict = {}
                    seq = ''
    
    #Convert list of entries to dataframe to easily save as file and use for analysis
    gbData_generated = pd.DataFrame(S_entries)
    
    return gbData_generated

def generate_processed_file(gbData_generated,OutputDirectory):

    gbData_generated.to_csv(OutputDirectory+"/sarsS_new.gb.processed.tsv", sep ='\t',index = False)

def parse_processed_file(GbProcessedFile):

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
        
        #Get first seq length that matches the longest sequence length identified
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

    print(combined_clusters)
    return cluster_centers,combined_clusters

def reassign_clusters(top_kmer_dict,cluster_centers,combined_clusters,gbData,OutputDirectory):

    for i in range(len(top_kmer_dict)):
        
        #If locus is not center locus for this cluster
        if(top_kmer_dict[i]['Locus'] != cluster_centers[top_kmer_dict[i]['Cluster']]):
            #Which locus does the cluster belong to
            center_locus = cluster_centers[top_kmer_dict[i]['Cluster']]
            #Reassign cluster value
            top_kmer_dict[i]['Cluster'] = combined_clusters[center_locus]
        else:
            #Reassign cluster value
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

def cluster_histogram(protein_clusters_df,OutputDirectory):

    x = list(protein_clusters_df['Cluster'])

    plt.hist(list(protein_clusters_df['Cluster']),bins = len(set(x)))
    axes = plt.gca()
    axes.set_ylim([0,1500])
    #Save figure to user specified output folder.
    plt.savefig(OutputDirectory+"/ClusterHistogram.png",bbox_inches='tight')
    plt.clf()


##Takes forever, probably won't do
##Create files for Cytoscape
def create_cluster_edge_list(combined_clusters,protein_clusters_df):

    clusters_edge_list = pd.DataFrame()

    for i in list(combined_clusters.values()):
        clusters_set = protein_clusters_df.loc[protein_clusters_df['Cluster'] == i]
        for j in range(len(clusters_set)-1):
            list_range = [item for item in range(j+1, len(clusters_set))]
            for k in list_range:
                clusters_edge_list = clusters_edge_list.append({'Source':clusters_set['Locus'].iloc[j],'Target':clusters_set['Locus'].iloc[k]},ignore_index=True)
    
    return clusters_edge_list


def output_cytoscape_files(clusters_edge_list,cluster_file,OutputDirectory):

    clusters_edge_list.to_csv(OutputDirectory+"/Clusters_EdgeList.csv", index = False)
    cluster_file.to_csv(OutputDirectory+"/Clusters_Metadata.noa", index = False)


def prepare_heatmap(protein_clusters_meta_df,OutputDirectory):

    #Combine all territories/region by country
    for i in range(len(protein_clusters_meta_df)):
        string = protein_clusters_meta_df['Location'][i].split(':')
        protein_clusters_meta_df['Location'][i] = string[0]

    #Output histogram of counts per country
    y = list(protein_clusters_meta_df['Location'])
    plt.hist(list(protein_clusters_meta_df['Location']),bins = len(set(y)))
    plt.xticks(rotation='vertical')
    plt.savefig(OutputDirectory+"/Location_Frequency.png")
    plt.clf()

    #Get frequency of just each country to normalize
    total_loc_df = protein_clusters_meta_df.groupby(['Location']).size().reset_index(name='TotalFrequency')

    #Get frequency of each occurance of location/cluster combo
    protein_clusters_loc_df = protein_clusters_meta_df.groupby(['Location', 'Cluster']).size().reset_index(name='ClusterFrequency')

    #Merge the 2 and calculate the % count per each location
    protein_clusters_loc_df = protein_clusters_loc_df.merge(total_loc_df,how='left', on='Location')
    protein_clusters_loc_df['NormFrequency'] = protein_clusters_loc_df['ClusterFrequency']/protein_clusters_loc_df['TotalFrequency']

    protein_clusters_loc_heatmap = pd.pivot_table(protein_clusters_loc_df,values='NormFrequency',index=['Location'], columns='Cluster')

    return protein_clusters_loc_heatmap

def generate_heatmap(protein_clusters_loc_heatmap,OutputDirectory):

    #Generate heatmap
    sns_plot = sns.heatmap(protein_clusters_loc_heatmap, cmap="YlGnBu")
    fig = sns_plot.get_figure()
    fig.savefig(OutputDirectory+"/ClusterLocation_Heatmap.png")

def main():

    #Generate argument parser and define arguments
    parser = defineArguments()
    args = parser.parse_args()
    
    GbOriginalFile = args.GbOriginalFile
    GbProcessedFile = args.GbProcessedFile
    OutputProcessedFile = args.OutputProcessedFile
    KmerSize = args.KmerSize
    NumKmersPerSequence = args.NumKmersPerSequence
    NumClusters = args.NumClusters
    OutputDirectory = args.OutputDirectory

    #Command expected to execute the script
    cmd = "--gb-original-file %s --output-processed-file %s --gb-processed-file %s --kmer-size %i --num-kmers-per-sequence %i --num-clusters %i --output-directory %s" % (GbOriginalFile,GbProcessedFile,OutputProcessedFile,KmerSize,NumKmersPerSequence,NumClusters,OutputDirectory)

    starttime = datetime.now()
    print("Sars_Protein_Clustering Start Time: ",starttime)

    #Ensure conflicting options are not provided
    if args.GbOriginalFile and args.GbProcessedFile:
        print("Conflicting options specified. Only specify original GenBank file or processed GenBank file to be analyzed. Exiting script.")
        sys.exit()

    if args.OutputProcessedFile and args.GbProcessedFile:
        print("Conflicting options specified. Processed GenBank file can only be generated when original GenBank file is provided without existing processed GenBank file. Exiting script.")
        sys.exit()

    #Parse and use original genBank file for analysis if provided
    if args.GbOriginalFile:
        gbData = parse_original_file(GbOriginalFile)

        #Output resulting processed file if specified
        if args.OutputProcessedFile:
            generate_processed_file(gbData,OutputDirectory)

    #Use processed genBank file for analysis if provided
    else:
        gbData = parse_processed_file(GbProcessedFile)

    #Prepare sequences
    AAhash = declare_AA_alphabet()
    gbData = replace_aas(AAhash,gbData)

    #Compute hash values for each kmer
    top_kmer_dict = compute_kmer_hashes(gbData,KmerSize,NumKmersPerSequence)

    print("Hashes computed: ",datetime.now())

    #Assign Clusters
    unique_hash = assign_clusters(top_kmer_dict)
    cluster_centers,combined_clusters = find_centers(top_kmer_dict)
    protein_clusters_df,protein_clusters_meta_df = reassign_clusters(top_kmer_dict,cluster_centers,combined_clusters,gbData,OutputDirectory)

    print("Clusters assigned: ",datetime.now())

    #Outputs
    cluster_histogram(protein_clusters_df,OutputDirectory)

    #Cytoscape: Takes a long time
    #cluster_edge_list = create_cluster_edge_list(combined_clusters,protein_clusters_df)
    #output_cytoscape_files(clusters_edge_list,cluster_file,OutputDirectory)

    #Generate heatmap
    protein_clusters_loc_heatmap = prepare_heatmap(protein_clusters_meta_df,OutputDirectory)
    generate_heatmap(protein_clusters_loc_heatmap,OutputDirectory)

    endtime = datetime.now()
    tdelta = endtime-starttime
    print("Total computation time: ",tdelta)

if __name__ == '__main__':
    main()
