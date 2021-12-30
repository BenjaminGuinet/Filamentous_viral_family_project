# Filamentous viral family project
------------
January 2022 
-------------

This project is a collaboration between the IMPI team in Tours (France) and the GEI team in Lyon (France). 
The objective of this project is to caracterize a new family of filamentous viruses.


We have in our possession 9 filamentous virus genomes, three of which are circularised. 
The idea of this pipeline is to :
- (i) characterize the genes of these viral genomes
- (ii) annotate them
- (iii) reconstruct the phylogeny of this group within the dsDNA virus clade

## Search for ORFs in the genomes 

Only open reading frames (ORFs) starting with a methionine and ending with a stop codon, with at least 50 amino acids and with minimal overlap (<23 nucleotides)
were considered as valuable candidates for being true ORFs.

**Programs** : 
- getorf (version EMBOSS:6.6.0.0)
- orf_filter.R 

**Output files**
```
- Virusname_prediction_option_1.pdf (Genomic plot of the filtred predicted ORFs for each contigs)
- Virusname_prediction_option_1.gff (Gff file with coordinates of each filtred ORFs)
- Virusname_prediction_option_1.fa (Fasta file with predicted filtred ORFs in AminoAcide format) 
- Virusname_getorf_option_1.fa (Fasta file with predicted ORFs from getorf)
```
------------------

## Gather and Cluster all homologous ORFs

- Knew LbFV-like dsDNA viruses : ```["LhFV_free","DFV_free","EfFV_free","PoFV_free","PcFV_free"]```
- Known dsDNA viruses : ```["AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","HzNV-1","HzNV-2","GbNV","OrNV","ToNV","DiNV","DmNV_kal","DmNV_tom","DmNV_esp","DmNV_mau","PmNV","HgNV","DhNV","GpSGHV","MdSGHV","LbFV","AmFV"]```

* Snakemake file : **Snakemake_Clustering**

* Script used : **Mmmseqs cluster, mmseqs createcsv,Mseqs2_clust_to_tab.py,Create_cluster_files.py**

Important file created : **ALL_Predicted_and_known_ORFs_cluster.tab** (contains all clusters of homologous loci 

**Clustering options :** 

- --cluster-mode 1 : Connected component clustering that uses the transitive connection to cover more distant homologs.
- --cov-mode 0 : Bidirectional coverage, where only sequences with overlapping sequence lengths greater than 30% of the longer of the two sequences are clustered, 
#(in this case always the viral sequence since the loci are defined by viral hits). 
- --evalue 0.0001 : To eliminate false positives. 
- --cluster-reassign : During the cascade clustering of Mmseqs2, as the representative of a cluster can change at each iteration, it can happen that some members that were already close to a cluster do not fulfill the clustering criteria anymore. We therefore correct this by reassigning these sequences.

- **Total number of Clusters : 1841**
- **Total number of Clusters (orphelins of known viruses removed) : 231/1841**
- **Total number of orphelins Clusters with new viruses ORFs : 142/231**


#All versus all 

/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search ALL_Predicted_and_known_ORFs_db ALL_Predicted_and_known_ORFs_db All_versus_All_Predicted_and_known_ORFs tpm_All_versus_All_Predicted_and_known_ORFs -s 7.5 -e 0.001 --threads 20

/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen' ALL_Predicted_and_known_ORFs_db ALL_Predicted_and_known_ORFs_db All_versus_All_Predicted_and_known_ORFs All_versus_All_Predicted_and_known_ORFs.m8

```
import networkx as nx
import pandas as pd 
import numpy as np

df=pd.read_csv("All_versus_All_Predicted_and_known_ORFs.m8",sep="\t")

df=df.loc[df['bits'].ge(50)]
# Create the graph from the dataframe
g = nx.Graph()

g.add_edges_from(df[['query','target']].itertuples(index=False))
new = list(nx.connected_components(g))

mapped =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new) for node in component}

cluster = pd.DataFrame({'Cluster': mapped.values(), 'Names':mapped.keys()})

cluster.to_csv("All_versus_All_Predicted_and_known_ORFs_clustered.m8",sep=";",index=False)

#Create binary table 

List_viruses_name=["AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","HzNV-1","HzNV-2","GbNV","OrNV","ToNV","DiNV","DmNV_kal","DmNV_tom","DmNV_esp","DmNV_mau","PmNV","HgNV","DhNV","GpSGHV","MdSGHV","LbFV","AmFV","PoFV","LhFV","PcFV","EfFV","DFV","CcFV1","CcFV2"]

#cluster="ALL_Predicted_and_known_ORFs_cluster.tab"
#cluster=pd.read_csv(cluster,header=0,sep="\t")
cluster['Names2'] =  cluster['Names'].str.replace('.*_\\+_|.*_-_','')

cluster['Names2'] = cluster['Names'].str.extract(f'({"|".join(List_viruses_name)})', expand=False)

#Load LbFV ORFs number table 

ORF_table = pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Clustering/LbFV_protein_ORF_names.txt",sep="\t")

cluster= cluster.merge(ORF_table,on='Names',how='outer')

cluster['Names2']=cluster['Names2'].replace('LbFV', np.nan)

cluster.Names2.fillna(cluster.New, inplace=True)
del cluster['New']
del cluster['Names']

List_ORFs=list(ORF_table['New'])

m = cluster['Names2'].isin(List_ORFs)

f = lambda x: [y for y in x if y]

cluster['ORFs_values'] = cluster['Cluster'].map(df['Names2'].where(m, '').groupby(cluster['Cluster']).agg(f))

cluster = cluster.explode('ORFs_values')

cluster = (pd.get_dummies(cluster.dropna(subset=['ORFs_values']), 
                     columns=['Names2'], 
                     prefix='', 
                     prefix_sep='')
        .groupby(['ORFs_values','Cluster'], as_index=False)
        .max())
        
cluster = cluster.loc[:, ~cluster.columns.str.startswith('LbFV_')]



cluster['ORF_numbers'] = cluster['ORFs_values'].str.replace('.*_ORF','').astype(int)
cluster.sort_values(by='ORF_numbers', ascending=True,inplace=True)

cluster=cluster[['ORFs_values', 'Cluster', 'LhFV', 'DFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','CpV', 'CuniNPV', 'DhNV', 'DiNV', 'DmNV_esp', 'DmNV_kal', 'DmNV_mau', 'DmNV_tom', 'GbNV', 'GpSGHV', 'HgNV', 'HzNV-1', 'HzNV-2', 'LdMNPV','MdSGHV', 'NeseNPV', 'OrNV',  'PmNV',  'ToNV','AcMNPV','AmFV']]

cluster[['ORFs_values', 'Cluster', 'LhFV', 'DFV','EfFV','PoFV','PcFV','CcFV1','CcFV2']]





--------------------

## Annotate proteins 

- Interproscan version : 5.53-87.0
- Database used : 
```[CDD-3.18,Coils-2.2.1,Gene3D-4.3.0,Hamap-2020_05,MobiDBLite-2.0,PANTHER-15.0,Pfam-34.0,PIRSF-3.10,PIRSR-2021_02,PRINTS-42.0,ProSitePatterns-2021_01,ProSiteProfiles-2021_01,SFLD-4,SMART-7.1,SUPERFAMILY-1.75,TIGRFAM-15.0]```

#Note : see this page if there there is issue between hmmer soft version and the database used : https://github.com/ebi-pf-team/interproscan/issues/173

Error example :

```
Error: File format problem in trying to open HMM file data/gene3d/4.3.0/gene3d_main.hmm.
Opened data/gene3d/4.3.0/gene3d_main.hmm.h3m, a pressed HMM file; but format of its .h3i file unrecognized
```

#Run snakemake on cluster 
```
nohup snakemake -j 8000  -s Snakemake_Interproscan  --cluster "sbatch -J {params.name} --mem {params.mem} -p normal -N 1 --cpus-per-task  {params.threads}  -o {params.out} -e {params.err}  " &> nohup_Interproscan_snakemake.out &

python3 /beegfs/home/bguinet/these_scripts_2/Add_interproscan_and_ontology.py -I ALL_Predicted_and_known_ORFs_intrproscan.tsv -c ALL_Predicted_and_known_ORFs_cluster.tab -o ALL_Predicted_and_known_ORFs_cluster_ontology.tab


```

* Snakemake file : **Snakemake_Interproscan**

* Script used : **Interproscan, Add_interproscan_and_ontology.py**

Important file created : **ALL_Predicted_and_known_ORFs_interproscan.tsv** (contains all ORFs with interproscan analysis)

--------------------

## Align & create phylogeny of clusters 

```
nohup snakemake -j 8000  -s Snakemake_Alignment_Phylogeny  --cluster "sbatch -J {params.name} --mem {params.mem} -p normal -N 1 --cpus-per-task  {params.threads}  -o {params.out} -e {params.err}  " &> nohup_Alignment_Phylogeny_snakemake.out &
```

* Snakemake file : **Snakemake_Alignment_Phylogeny**

* Script used : **macse_v2.05.jar, Iqtree-2.1.2**

Important file created :
- **{cluster_number}_AA.dna** (contains the codon alignment of the ORFs in each clusters )
- **{cluster_number}_AA.dna.treefile** (contains the phylogenetic tree of the ORFs in each clusters > 2 ORFs )


The FLV clusters ontologies 


tab=pd.read_csv("ALL_Predicted_and_known_ORFs_cluster_ontology.tab",sep=";")

tab=tab.loc[tab['Names'].str.contains("FV")]
tab=tab.loc[~tab['Names'].str.contains("AmFV")]


ok

tab[['Clustername','Gene names','Protein names','PANTHER_Interpro_description','CDD_Interpro_description','Pfam_Interpro_description','TIGRFAM_Interpro_description','ProSiteProfiles_Interpro_description','Pfam_Signature_description','PIRSF_Interpro_description','Pfam_Interpro_description','SMART_Signature_description]]

tab = tab[tab['Gene names'].notna() |tab['Protein names'].notna() | tab['PANTHER_Interpro_description'].notna() |tab['Hamap_Interpro_description'].notna() |tab['CDD_Interpro_description'].notna() |tab['Pfam_Interpro_description'].notna() |tab['TIGRFAM_Interpro_description'].notna() |tab['ProSiteProfiles_Interpro_description'].notna() | tab['Pfam_Signature_description'].notna() | tab['PIRSF_Interpro_description'].notna()|tab['Pfam_Interpro_description'].notna() | tab['SMART_Signature_description'].notna()  ]

Hamap_Interpro_description



tab[['Gene names','Protein names','PANTHER_Interpro_description','Hamap_Interpro_description']
#Keep only cluster with at leats one ontology assignation


tab.loc[Cluster','Names','Gene names', 'Protein names','PANTHER_Interpro_description','Hamap_Interpro_description','CDD_Interpro_description','Pfam_Interpro_description','TIGRFAM_Interpro_description']]

_________________


## Concatenate phylogeny analysis 

#First we need to create new cluster alignment file with same tips names 

```
import os 
import re 
from Bio import SeqIO 

viruses = ["AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","HzNV-1","HzNV-2","GbNV","OrNV","ToNV","DiNV","DmNV_kal","DmNV_tom","DmNV_esp","DmNV_mau","PmNV","HgNV","DhNV","GpSGHV","MdSGHV","LbFV","AmFV","PoFV","LhFV","PcFV","EfFV","DFV","CcFV1","CcFV2"]

directory="/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/"
for filename in os.listdir(directory):
  if "_AA.dna" in filename:
    clustername= re.sub("_AA.dna","",filename)
    nb_viruses=0
    record=record_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    nb_viruses=len(record)
    if nb_viruses > 4:
      with open(directory+clustername+"_AA2.dna","w") as output:
          viruses_added_to_cluster=[]
          for i in record:
            for viruse in viruses:
              if viruse in record[i].id:
                if viruse not in viruses_added_to_cluster:
                  viruses_added_to_cluster.append(viruse)
                  print('>',viruse,sep="",file=output)
                  print(record[i].seq,file=output)
```

Then we need to concatenate those files :

```
cd /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/
perl /beegfs/home/bguinet/these_scripts_2/catfasta2phyml.pl -f --concatenate --verbose *_AA2.dna  > /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/Concatenated_sequences.aln  2> partitions.txt

#Run the phylogeny

/beegfs/data/bguinet/Bguinet_conda/bin/iqtree -s /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/Concatenated_sequences.aln -spp /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/partitions.tab -m TEST -alrt 1000  -bb 1000  -nt 25


```

/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/*_AA.dna 


