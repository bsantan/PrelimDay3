#######################################################################
##
##
## 
######################################################################
#
#The Sars_Protein_Clustering.py is an algorithm that takes a GenBank file consisting of complete and partial sequences of viral isolates and clusters the surface glycoprotein ("S") sequences according to sequence similarity. The program can either be input a complete GenBank file with all records of a given query, after which it will extract the "S" proteins of interest, or a processed file with the "S" protein sequences already extracted. The algorithm will apply a Kmer based method to assign each variant to a cluster, and will output the result of these cluster assignments to a csv file. The program will also expand the results of the analysis with several visualization techniques across time and space of the protein variants given.

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
import sys
import dateutil.parser

#Libraries for validation only
import difflib
from difflib import SequenceMatcher

#Define arguments for each required and optional input
def defineArguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ##Need to add logic for how to handle these entries
    #When provided, will use output as gbData
    parser.add_argument("--gb-original-file",  dest="GbOriginalFile",required=False,help="GbOriginalFile")
    parser.add_argument("--output-processed-file",  dest="OutputProcessedFile",required=False,action="store_true",help="OutputProcessedFile")
    #When provided, will not perform parsing on original genBank file
    parser.add_argument("--gb-processed-file",  dest="GbProcessedFile",required=False,help="GbProcessedFile")
    #For validation purposes
    parser.add_argument("--og-protein-fasta",  dest="OgProteinFasta",required=False,help="OgProteinFasta")
    #Parameters for algorithm
    parser.add_argument("--kmer-size", dest="KmerSize",required=False,default=10,help="KmerSize",type=int)
    parser.add_argument("--num-kmers-per-sequence", dest="NumKmersPerSequence",required=False,default=10,help="NumKmersPerSequence",type=int)
    parser.add_argument("--output-directory",  dest="OutputDirectory",required=True,help="OutputDirectory")

    return parser

#Parses complete GenBank file to result in LocusID, location, collectionDate, and AAsequence of only S proteins. Resulting df will be used throughout clustering algorithm. 
def parse_original_file(GbOriginalFile):

    f = open(GbOriginalFile, "r")

    S_entries = []
    #Initialize the dict once
    entry_dict = {}
    desired_gene = False
    getting_protein = False
    seq = ''

    for i,line in enumerate(f):
        #Manage empty lines
        line = line.strip()
        #Extract tab delimited values
        values = line.split()
        if values:
            #Get all necessary items for this entry
            if values[0].lower() == 'accession':
                entry_dict['LOCUS'] = values[1]
            if 'country' in values[0].lower():
                #For USA locations with state listed, or other country where region comes after colon
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
            
            #Only search for "translation" string if gene="s" was found
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

#Outputs result of parse_original_file to a csv file.
def generate_processed_file(gbData_generated,OutputDirectory):

    gbData_generated.to_csv(OutputDirectory+"/sarsS_new.gb.processed.tsv", sep ='\t',index = False)

#Parses a processed GenBank file that already contains proteins of interest only. Resulting df will be used iin clustering algorithm. 
def parse_processed_file(GbProcessedFile):

    gbData_og = pd.read_csv(GbProcessedFile, sep = '\t')

    return gbData_og

#Parses original "reference" protein sequence from a fasta file to compare results to.
def parse_protein_file(OgProteinFasta):

    read = ''
    f = open(OgProteinFasta)
    lines = f.readlines()
    for i in range(len(lines)):
        if not lines[i].startswith(">"):
            read += lines[i].strip()

    return read

#Defines how amino acid characters will be converted for kmer comparison.
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

#Replaces amino acid characters among all sequences, resulting in gbData which will be used in clustering algorithm.
def replace_aas(AAhash,gbData_og):

    gbData = copy.deepcopy(gbData_og)

    for i in range(len(gbData_og)):
        for j in (AAhash.keys()):
            gbData['AAsequence'][i] = gbData_og['AAsequence'][i].replace(j,AAhash[j])

    return gbData

#Computes a hash value for each kmer and outputs the hash along with kmer sequence, Locus ID, index, and length to a df.
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
            try:
                #Gets ASCII value for each character
                kmer_hash += ord(kmer[j])**(KmerSize - (j+1))
            except: 
                #Handle case where kmer size specified is longer than AA sequence
                print("Kmer size specified is longer than AAsequence of length: "+str(len(AAsequence)))
                sys.exit()
        #Gets ASCII value for each character
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
    
    top_kmer_dict = pd.DataFrame(top_kmer_dict)

    return top_kmer_dict

#Assigns all unique kmer hash values to a cluster by enumerating them, and combining those that are the same.
def assign_clusters(top_kmer_dict):

    #List that will keep track of already identified hash values vs. new ones
    unique_hash = []
    top_kmer_dict['Cluster'] = None

    for i in range(len(top_kmer_dict)):
        #Add new hash to list if not already there, and assign cluster
        if(top_kmer_dict['Kmer_hash'][i] not in unique_hash):
            unique_hash.append(top_kmer_dict['Kmer_hash'][i])
            top_kmer_dict['Cluster'][i] = (unique_hash.index(top_kmer_dict['Kmer_hash'][i])+1)
        #Only assign cluster if hash is already in list
        else:
            top_kmer_dict['Cluster'][i] = (unique_hash.index(top_kmer_dict['Kmer_hash'][i])+1)

    return top_kmer_dict

#Identifies center protein of each unique cluster/kmer hash value by getting the longest sequence, and assigns a new cluster value to that. Will also find the sequence similarity by comparing to og sequence or reference protein fasta.
def find_centers(top_kmer_dict,gbData_og,og_sequence,OutputDirectory):

    #Get all unique clusters assigned
    unique_clusters = top_kmer_dict.Cluster.unique()
    cluster_centers = {}

    for i in unique_clusters:
        #Subset df by only this cluster
        cluster_kmer_df = top_kmer_dict.loc[top_kmer_dict['Cluster'] == i].reset_index()
        #Identify the longest protein sequence from that cluster
        top_length = max(cluster_kmer_df['Seq_length'])
        
        #Get first seq length that matches the longest sequence length identified
        for j in range(len(cluster_kmer_df)):
            if(cluster_kmer_df['Seq_length'][j] == top_length):
                cluster_centers[cluster_kmer_df['Cluster'][j]] = cluster_kmer_df['Locus'][j]
                #Only use first sequence of that length
                break

    #Get all unique centers to assign new cluster to 
    all_centers = pd.DataFrame(cluster_centers.values())
    all_centers = all_centers.drop_duplicates()
    all_centers = all_centers.reset_index(drop=True)

    #Dictionary to represent the new assigned cluster for each center [Locus, Cluster]
    combined_clusters = {}
    #Dictionary to represent the similarity with the original S protein for each center of a Cluster [Cluster, Similarity]
    center_similarities = {}

    clust = 1
    for i in all_centers[0]:
        #Assign new cluster to each unique center
        combined_clusters[i] = clust
        #Evaluate sequence similarity for validation purposes while looping through all centers
        sequence = gbData_og.loc[gbData_og['LOCUS'] == i,'AAsequence'].iloc[0]
        similarity = difflib.SequenceMatcher(None, og_sequence, sequence).ratio()
        center_similarities[clust] = similarity
        clust += 1

    #Convert similarities to dataframe
    center_similarities = pd.DataFrame(list(center_similarities.items()))
    center_similarities.columns = ['Cluster', 'Similarity']

    #Write sequence similarities to table for analysis
    center_similarities.to_csv(OutputDirectory+"/Clusters_Sequence_Similarities.csv", index = False)

    return cluster_centers,combined_clusters,center_similarities

#Uses new cluster values to reassign each locus to a cluster.
def reassign_clusters(top_kmer_dict,cluster_centers,combined_clusters,gbData,OutputDirectory):

    #Get all unique loci from df
    unique_locus = top_kmer_dict.Locus.unique()
    #Will contain each unique locus and its new cluster assignment 
    protein_clusters_df = pd.DataFrame()

    #For each locus, identify most common cluster assignment among kmers and assign the locus to that
    for i in unique_locus:
        #Extract all entries from this locus
        locus_set = top_kmer_dict.loc[top_kmer_dict['Locus'] == i]
        #Use the original cluster with more kmers assigned to identify the center if there are more than 1 for this locus
        if(len(locus_set.Cluster.unique()) > 1):
            selected_cluster = locus_set['Cluster'].value_counts().idxmax()
        else: selected_cluster = locus_set['Cluster'].iloc[0]
        #Reassign cluster
        protein_clusters_df = protein_clusters_df.append({'Locus':i,'Seq_length':locus_set['Seq_length'].iloc[0],'Cluster':selected_cluster},ignore_index=True)
    
    unique_clusters = protein_clusters_df.Cluster.unique()

    #Assign new cluster labels that are continuos
    #Set counter for new cluster label
    new_cluster = 1
    for i in unique_clusters:
        protein_clusters_df.loc[protein_clusters_df.Cluster == i, 'Cluster'] = new_cluster
        new_cluster += 1

    #Convert Cluster labels to integers for readability
    protein_clusters_df['Cluster'] = protein_clusters_df['Cluster'].astype(int)

    #Add metadata
    protein_clusters_meta_df = pd.concat([protein_clusters_df,gbData[['Location', 'collectionDate']].reset_index(drop=True)],axis=1)

    protein_clusters_df[['Locus','Cluster']].to_csv(OutputDirectory+"/Clusters.csv", index = False)

    return protein_clusters_df,protein_clusters_meta_df

#Creates histogram of # of loci assigned to each cluster.
def generate_cluster_histogram(protein_clusters_df,OutputDirectory):

    #Get unique set of clusters
    x = list(protein_clusters_df['Cluster'])

    plt.hist(list(protein_clusters_df['Cluster']),bins = len(set(x)))
    axes = plt.gca()
    #axes.set_ylim([0,1500])
    plt.xlabel("Cluster")
    plt.ylabel("# Variants Assigned")
    plt.title('Number of Protein Variants per Cluster', fontsize=10)
    plt.xticks(rotation='vertical')
    plt.autoscale()
    #Save figure to user specified output folder.
    plt.savefig(OutputDirectory+"/ClusterHistogram.png",bbox_inches='tight')
    plt.clf()

#Creates histogram of number of samples existing across all dates.
def eval_location_frequency(protein_clusters_meta_df,OutputDirectory):

    #Initialize column for general location.
    protein_clusters_meta_df['General_Location'] = None

    #Combine all territories/region by country
    for i in range(len(protein_clusters_meta_df)):
        string = protein_clusters_meta_df['Location'][i].split(':')
        protein_clusters_meta_df['General_Location'][i] = string[0]

    #Output histogram of counts per country
    y = list(protein_clusters_meta_df['General_Location'])
    #Bin by total number of unique locations.
    plt.hist(list(protein_clusters_meta_df['General_Location']),bins = len(set(y)))
    plt.xticks(rotation='vertical')
    plt.xlabel("Location")
    plt.ylabel("# Samples")
    plt.title('Number of Samples per Location', fontsize=10)
    plt.autoscale()
    plt.savefig(OutputDirectory+"/Location_Frequency.png",bbox_inches="tight")
    plt.clf()

#Identifies average date for each cluster.
def get_avg_dates(protein_clusters_meta_df,OutputDirectory):

    #Convert all dates to datetime
    for i in range(len(protein_clusters_meta_df)):
        protein_clusters_meta_df['collectionDate'][i] = dateutil.parser.parse(protein_clusters_meta_df['collectionDate'][i])

    #Find average date per cluster
    #Get only unique clusters
    unique_final_clusters = protein_clusters_meta_df.Cluster.unique()
    #Create dict for avg date of each cluster
    avg_date = {}

    for i in unique_final_clusters:
        clusters_set = protein_clusters_meta_df.loc[protein_clusters_meta_df['Cluster'] == i]
        upper_d = max(clusters_set['collectionDate'])
        lower_d = min(clusters_set['collectionDate'])
        avg_d = lower_d+((upper_d-lower_d)//2)
        
        avg_date[i] = avg_d
        
    #Assign average dates to each cluster
    #Initialize column for average date
    protein_clusters_meta_df['avgCollectionDate'] = None

    for i in range(len(protein_clusters_meta_df)):
        protein_clusters_meta_df['avgCollectionDate'][i] = avg_date[protein_clusters_meta_df['Cluster'][i]]

    protein_clusters_meta_df.to_csv(OutputDirectory+"/Clusters_Meta.csv", index = False)

    return protein_clusters_meta_df

#Generates scatter plot of average date for each cluster.
def plot_clusters_avg_date(protein_clusters_meta_df,OutputDirectory):

    #Sort clusters by average date, convert to string to be categorical.
    protein_clusters_meta_df.sort_values(by=['avgCollectionDate'],inplace=True,ascending=True)
    protein_clusters_meta_df['Cluster'] = protein_clusters_meta_df['Cluster'].astype(str)

    #Plot of average date collected for each cluster of variants
    plt.scatter(protein_clusters_meta_df['Cluster'],protein_clusters_meta_df['avgCollectionDate'])
    plt.xlabel("Cluster")
    plt.ylabel("Average Date")
    plt.title('Average Date Per Cluster', fontsize=10)
    plt.xticks(rotation='vertical')

    plt.savefig(OutputDirectory+"/Clusters_Average_Date.png",bbox_inches="tight")
    plt.clf()

#Generates histogram of samples collected at each date.
def plot_samples_over_time(protein_clusters_meta_df,OutputDirectory):

    plt.hist(protein_clusters_meta_df['collectionDate'],bins=50)
    plt.xlabel("Date Collected")
    plt.ylabel("# Samples")
    plt.title('Number of Samples Collected Over Time', fontsize=10)
    plt.xticks(rotation='vertical')
    plt.autoscale()
    plt.savefig(OutputDirectory+"/Date_Frequency.png",bbox_inches="tight")
    plt.clf()

#Calculates fraction of samples at each location assigned to each cluster and stores as a table. A table for all locations and non-US locations only will be generated.
def prepare_heatmap(protein_clusters_meta_df,OutputDirectory):

    #Get frequency of just each country across all clusters to normalize (TotalFrequency)
    total_loc_df = protein_clusters_meta_df.groupby(['Location']).size().reset_index(name='TotalFrequency')

    #Get frequency of each country for each cluster (ClusterFrequency)
    protein_clusters_loc_df = protein_clusters_meta_df.groupby(['Location', 'Cluster']).size().reset_index(name='ClusterFrequency')

    #Merge the 2 for calculation purposes (NormFrequency), and include avgCollectionDate
    protein_clusters_loc_df = protein_clusters_loc_df.merge(total_loc_df,how='left', on='Location')
    protein_clusters_loc_df = protein_clusters_loc_df.merge(protein_clusters_meta_df[['Cluster','avgCollectionDate']],how='left', on='Cluster')
    #Calculate normalized frerquency for each location
    protein_clusters_loc_df['NormFrequency'] = protein_clusters_loc_df['ClusterFrequency']/protein_clusters_loc_df['TotalFrequency']

    #Sort clusters by average date
    protein_clusters_loc_df.sort_values(by=['avgCollectionDate'],inplace=True,ascending=True)
    #Convert Cluster name to string to use these as categorical values
    protein_clusters_loc_df['Cluster'] = protein_clusters_loc_df['Cluster'].astype(str)

    #Generate dataframe for all samples not in the US
    protein_clusters_notUS_df = protein_clusters_loc_df[~protein_clusters_loc_df.Location.str.contains("USA")]

    #Generate heatmap for normalized frequency of each county per cluster
    protein_clusters_loc_heatmap = pd.pivot_table(protein_clusters_loc_df,values='NormFrequency',index=['Location'], columns='Cluster')
    #Also generate table for locations not in the US
    protein_clusters_notUS_heatmap = pd.pivot_table(protein_clusters_notUS_df,values='NormFrequency',index=['Location'], columns='Cluster')

    #Reorder so that clusters go by ascending average date 
    all_clusters = protein_clusters_loc_df['Cluster'].drop_duplicates()
    all_clusters = list(all_clusters.reset_index(drop=True))
    protein_clusters_loc_heatmap = protein_clusters_loc_heatmap.reindex(all_clusters, axis=1)
    #Repeat for non US locations
    all_clusters_notUS = protein_clusters_notUS_df['Cluster'].drop_duplicates()
    all_clusters_notUS = list(all_clusters_notUS.reset_index(drop=True))
    protein_clusters_notUS_heatmap = protein_clusters_notUS_heatmap.reindex(all_clusters_notUS, axis=1)

    #Write heatmap pivot table to csv for analysis
    protein_clusters_loc_heatmap.to_csv(OutputDirectory+"/Clusters_Heatmap_All.csv")
    protein_clusters_notUS_heatmap.to_csv(OutputDirectory+"/Clusters_Heatmap_notUS.csv")

    return protein_clusters_loc_heatmap,protein_clusters_notUS_heatmap

#Generates a heatmap based on above function for fraction of samples in each location that belong to each cluster, ordered by avg date collected. A heatmap for all locations and non-US locations only will be generated.
def generate_heatmap(protein_clusters_loc_heatmap,loc,OutputDirectory):

    #Generate heatmap
    sns_plot = sns.heatmap(protein_clusters_loc_heatmap, cmap="flare")
    plt.title('Fraction of Variants Assigned to Clusters for Each Location')
    plt.autoscale()
    plt.yticks(rotation='horizontal')
    fig = sns_plot.get_figure()
    fig.savefig(OutputDirectory+"/ClusterLocation_Heatmap"+loc+".png",bbox_inches="tight")
    plt.clf()

#Uses calculated similarity of each center protein to output a plot of the top 10 largest clusters and their similarity. For validation purposes only. 
def output_similarity_metrics(center_similarities,og_sequence,protein_clusters_meta_df,OutputDirectory):

    #Find the top 10 largest clusters
    cluster_frequencies = protein_clusters_meta_df.groupby(['Cluster']).size().reset_index(name='TotalClusterFrequency')
    cluster_frequencies = cluster_frequencies.nlargest(10, ['TotalClusterFrequency'])

    #Sort the top 10 clusters by average collection date
    center_similarities = center_similarities[center_similarities['Cluster'].isin(list(cluster_frequencies['Cluster']))]
    center_similarities['Cluster'] = center_similarities['Cluster'].astype(str)
    center_similarities = center_similarities.merge(protein_clusters_meta_df[['Cluster','avgCollectionDate']],how='left', on='Cluster')
    center_similarities.sort_values(by=['avgCollectionDate'],inplace=True,ascending=True)

    #Plot sequence similarity
    plt.scatter(center_similarities['Cluster'],center_similarities['Similarity'])
    plt.xlabel("Cluster")
    plt.ylabel("Sequence Similarity")
    plt.title('Similarity of Cluster Center Proteins with Original Sequence', fontsize=10)
    plt.xticks(rotation='vertical')
    plt.autoscale()
    plt.savefig(OutputDirectory+"/Sequence_Similarity_top10.png",bbox_inches="tight")
    plt.clf()

def main():

    #Manage warning about assigning new values to each dataframe
    pd.options.mode.chained_assignment = None

    #Generate argument parser and define arguments
    parser = defineArguments()
    args = parser.parse_args()
    
    GbOriginalFile = args.GbOriginalFile
    GbProcessedFile = args.GbProcessedFile
    OutputProcessedFile = args.OutputProcessedFile
    OgProteinFasta = args.OgProteinFasta
    KmerSize = args.KmerSize
    NumKmersPerSequence = args.NumKmersPerSequence
    OutputDirectory = args.OutputDirectory

    #Command expected to execute the script
    cmd = "--gb-original-file %s --output-processed-file %s --gb-processed-file %s --og-protein-fasta %s --kmer-size %i --num-kmers-per-sequence %i --output-directory %s" % (GbOriginalFile,GbProcessedFile,OutputProcessedFile,OgProteinFasta,KmerSize,NumKmersPerSequence,OutputDirectory)

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
        gbData_og = parse_original_file(GbOriginalFile)

        #Output resulting processed file if specified
        if args.OutputProcessedFile:
            generate_processed_file(gbData_og,OutputDirectory)

    #Use processed genBank file for analysis if provided
    else:
        gbData_og = parse_processed_file(GbProcessedFile)

    #Read S protein reference file
    og_sequence = parse_protein_file(OgProteinFasta)

    #Prepare sequences
    AAhash = declare_AA_alphabet()
    gbData = replace_aas(AAhash,gbData_og)

    #Compute hash values for each kmer
    top_kmer_dict = compute_kmer_hashes(gbData,KmerSize,NumKmersPerSequence)

    print("Hashes computed: ",datetime.now())

    #Assign Clusters
    top_kmer_dict = assign_clusters(top_kmer_dict)
    cluster_centers,combined_clusters,center_similarities = find_centers(top_kmer_dict,gbData_og,og_sequence,OutputDirectory)
    protein_clusters_df,protein_clusters_meta_df = reassign_clusters(top_kmer_dict,cluster_centers,combined_clusters,gbData,OutputDirectory)

    print("Clusters assigned: ",datetime.now())

    #Output ClusterHistogram.png
    generate_cluster_histogram(protein_clusters_df,OutputDirectory)

    #Evaluate clusters across time
    eval_location_frequency(protein_clusters_meta_df,OutputDirectory)
    protein_clusters_meta_df = get_avg_dates(protein_clusters_meta_df,OutputDirectory)

    #Output Clusters_Average_Date.png, Date_Frequency.png
    plot_clusters_avg_date(protein_clusters_meta_df,OutputDirectory)
    plot_samples_over_time(protein_clusters_meta_df,OutputDirectory)

    #Generate heatmap of cluster counts per location in order of average date of each cluster
    protein_clusters_loc_heatmap,protein_clusters_notUS_heatmap = prepare_heatmap(protein_clusters_meta_df,OutputDirectory)

    #Output ClusterLocation_HeatmapAll.png, ClusterLocation_HeatmapNotUS.png
    generate_heatmap(protein_clusters_loc_heatmap,"All",OutputDirectory)
    generate_heatmap(protein_clusters_notUS_heatmap,"NotUS",OutputDirectory)

    #Output Sequence_Similarity_top10.png
    output_similarity_metrics(center_similarities,og_sequence,protein_clusters_meta_df,OutputDirectory)

    endtime = datetime.now()
    tdelta = endtime-starttime
    print("Total computation time: ",tdelta)

if __name__ == '__main__':
    main()
