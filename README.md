# Filamentous_viral_family_project

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

_________________


## Concatenate phylogeny analysis 

#First we need to create new cluster alignment file with same tips names 

```
import os 
from Bio import SeqIO 

for filename in os.listdir(directory):
  if "_AA.dna" in filename:
    clustername= re.sub("_AA.dna","",filename)
    nb_viruses=0
    record=record_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    nb_viruses=len(record)
    if nb_viruses > 4:
      with open("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/"+clustername+"_AA2.dna","w") as output:
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
/beegfs/data/bguinet/Bguinet_conda/bin/iqtree -s /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/Concatenated_sequences.aln -spp /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/partitions.tab -m TEST -alrt 1000  -bb 1000  -nt 20
```

/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_alignment/*_AA.dna 


