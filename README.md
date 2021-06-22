# Clustering Sars-CoV-2 Protein Sequences

The Sars_Protein_Clustering.py is an algorithm that takes a GenBank file consisting of complete and partial sequences of viral isolates and clusters the surface glycoprotein ("S") sequences according to sequence similarity. The program can either be input a complete GenBank file with all records of a given query, after which it will extract the "S" proteins of interest, or a processed file with the "S" protein sequences already extracted. When a complete GenBank file is given, the program will only identify S protein sequences, however with an already processed file the algorithm is not protein specific.

## Getting Started

These instructions will provide the necessary environments, programs with installation instructions, and input files in order to run the Sars_Protein_Clustering.py script.

### Prerequisites
The following software packages must be installed:
```
- Python (any version later than 2.7)
```

## Running the Script

To run the Sars_Protein_Clustering.py script with a complete GenBank file, execute the following command, where GBORIGINALFILE is a complete GenBank file from an NCBI query, OGPROTEINFASTA is the reference S protein fasta file, and OUTPUTDIRECTORY is the desired output directory:

```
Sars_Protein_Clustering.py --gb-original-file GBORIGINALFILE --og-protein-fasta OGPROTEINFASTA --output-directory OUTPUTDIRECTORY
```

To run the script with a processed GenBank file including only protein sequences of interest, execute the following command, where GBPROCESSEDFILE is a tab-delimited file with the columns LOCUS, Location, collectionDate, and AAsequence:

```
Sars_Protein_Clustering.py --gb-processed-file GBPROCESSEDFILE --og-protein-fasta OGPROTEINFASTA --output-directory OUTPUTDIRECTORY
```

*Note that the specified output directory must already exist

The commands --gb-original-file and --gb-processed-file cannot be used simultaneously. 


### Required Input Files

To run the Sars_Protein_Clustering.py script, a complete GenBank file OR a processed GenBank file are required, as well as a FASTA file containing the reference "S" protein sequence.

The format of the complete GenBank file is described in detail here: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html 

The format of the processed GenBank file is:
```
LOCUS	Location	collectionDate	AAsequence
MT079851	China	2020-01-22	MFVFLVLLPLVSSQCVN
```

The format of the FASTA file of the reference proteins sequence is:
```
>FASTA_ID
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHV
SGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPF
```

## Expected Outputs

The Sars_Protein_Clustering.py script will generate the following file when the --output-processed-file is enabled (default). 

### sarsS_new.gb.processed.tsv
This file mimics the structure of the processed GenBank (GBPROCESSEDFILE) above.

To disable this option, include the following command: 

```
--output-processed-file false
```

The Sars_Protein_Clustering.py script will always generate the following files:

### Clusters.csv
A table containing each locus and the corresponding cluster assignment in the following format:
```
Locus	Cluster	
MT079851	1	
```

The program will also output a histogram of this information as ClusterHistogram.png. 

### Clusters_Meta.csv
A table containing each locus and the corresponding cluster assignment along with metadata in the following format:	
```
Cluster	Locus	Seq_length	Location	collectionDate	General_Location	avgCollectionDate
1	MT079851	1273	China	1/22/20 0:00	China	7/22/20 12:00
```

### Clusters_Heatmap_All.csv and Clusters_Heatmap_notUS.csv
A table containing the normalized frequency of cluster assignments (top row) for each location (left column), where clusters are sorted by average date collected among samples assigned. The file with "_All" includes all locations, and the table with "_notUS" includes only locations excluding "USA" in the designated collection location. Both tables are in the following format:
```
Location	4	2	5	8	6	9	1	7	3
China							1.00		
```						

The program will also output heatmaps of this information for visualization as ClusterLocation_HeatmapAll.png and ClusterLocation_HeatmapNotUS.png.

### Clusters_Sequence_Similarities.csv
A table showing the sequence similarity of each cluster center with the reference protein sequence provided as OGPROTEINFASTA in the following format:
```
Cluster	Similarity
1	1
```

The program will also output a plot of this information including only the largest 10 clusters in order of average date collected as Sequence_Similarity_top10.png. 

### Other plots for visualization

```
- Clusters_Average_Date.png: a plot showing the average date collected for samples assigned to each cluster 
- Date_Frequency.png: a histogram showing the number of samples collected across each date
- Location_Frequency.png: a histogram showing the number of samples collected across each general location (most often country)
```


### Optional Arguments

By default, the program will use kmers of length 10 to do pairwise comparisons between protein sequences from different samples. However you can change the kmer length by including the following optional argument. If the size specified is longer than the shortest identified protein sequence, the script will exit:

```
--kmer-size KMERSIZE (default 10)
```

Additionally, the program will select the 10 kmers from each sequence to do the pairwise comparisons. However you can change the number of kmers selected by including the following optional argument:

```
--num-kmers-per-sequence NUMKMERSPERSEQUENCE (default 10)
```
