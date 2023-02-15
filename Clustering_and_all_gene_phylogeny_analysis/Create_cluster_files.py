import sys 
from Bio import SeqIO
from Bio import Entrez
import pandas as pd 
from Bio import SeqIO
import os
import argparse


# Print out a message when the program is initiated. traitrise
print('----------------------------------------------------------------------------------------------\n')
print('                        Assign ORFs sequences to cluster created previously           .\n')
print('----------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Assign ORFs sequences to cluster files')
parser.add_argument("-c", "--cluster", help="The cluster file in tabular separator format")
parser.add_argument("-o", "--Outdir", help="The output directory of the clusters")
parser.add_argument("-AA", "--AA", help="The file with AA ORFs")
parser.add_argument("-DNA", "--DNA", help="The file with DNA ORFs")
args = parser.parse_args()


#Usage example python3 /beegfs/home/bguinet/these_scripts_2/Create_cluster_files.py -c /beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab -AA /beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs.fa -o /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs/ -DNA /beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_known_ORFs.fna
def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

cluster_file=args.cluster
AA_virus_ORFs=args.AA
DNA_virus_ORFs=args.DNA
Output_directory=args.Outdir


record_dict = to_dict_remove_dups(SeqIO.parse(AA_virus_ORFs, "fasta"))
cluster_tab=pd.read_csv(cluster_file,sep=";")
print(cluster_tab)
cluster_tab = cluster_tab.drop_duplicates(subset = 'Names', keep = 'first')

#Write Cluster files in Amino Acide format
for _, cluster in cluster_tab.groupby('Cluster'):
	print(str(cluster['Cluster'].unique()[0]))
	for row1, row2  in zip(cluster['Cluster'], cluster['Names']):
		filename=Output_directory+row1+".aa"
		with open(filename, "a+") as output_handle:
				print(">",record_dict[row2].id, file=output_handle,sep="")
				print(record_dict[row2].seq,file=output_handle)
	print(row1, "Viral ORFs recovery done")

#Write Cluster files in Nucleotides format

#Add the new viruses sequences :

list_new_viruses=['PoFV','PcFV','EfFV','DFV','LhFV','CcFV1','CcFV2']

for viruses in list_new_viruses:
	sub_cluster= cluster_tab.loc[cluster_tab['Names'].str.contains(viruses)]
	sub_cluster['Names2'] = sub_cluster['Names'].str.replace('-','_')
	sub_cluster['Names2'] = sub_cluster['Names2'].str.replace('+','_')
	sub_cluster['Names2'] = sub_cluster['Names2'].str.replace('___.*','')
	sub_cluster[['Scaffold', 'start', 'end']] = sub_cluster['Names2'].str.rsplit('_', n=2, expand=True)
	sub_cluster['start']=sub_cluster['start'].astype(int)
	sub_cluster['end']=sub_cluster['end'].astype(int)
	for _, sub_sub_cluster in sub_cluster.groupby('Cluster'):
		if viruses == "PoFV":
			name= "Porseoliae"
		elif viruses == "PcFV":
			name = "Pconcolor"
		elif viruses == "EfFV":
			name = "Eformosa"
		elif viruses == "DFV":
			name = "DFV"
		elif viruses == "LhFV":
			name= "LhFV"
		elif viruses == "CcFV1":
			name= "CcFV1"
		elif viruses == "CcFV2":
			name= "CcFV2"
		record_dict = SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+viruses+"_free/"+name+"_LbFV_free.fna", "fasta"))
		filename=Output_directory+sub_sub_cluster['Cluster'].iloc[0]+".dna"
		with open(filename, "a") as output:
			for index, row in sub_sub_cluster.iterrows():
				if "+" in row['Names']:
					seq=str(record_dict[row['Scaffold']][row['start']-1: row['end']].seq)
					print(">",row['Names'],sep="",file=output)
					print(seq,file=output)
				elif "-" in row['Names']:
					seq=str(record_dict[row['Scaffold']][row['start']-1: row['end']].seq.reverse_complement())
					print(">",row['Names'],sep="",file=output)
					print(seq,file=output)




List_viruses_name=["AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","HzNV-1","HzNV-2","GbNV","OrNV","ToNV","DiNV","DmNV_kal","DmNV_tom","DmNV_esp","DmNV_mau","PmNV","HgNV","DhNV","GpSGHV","MdSGHV","LbFV","AmFV","PoFV","LhFV","PcFV","EfFV","DFV","CcFV1","CcFV2","WSSV","CoBV"]
List_viruses_name_without_outgroupp=["AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","HzNV-1","HzNV-2","GbNV","OrNV","ToNV","DiNV","DmNV_kal","DmNV_tom","DmNV_esp","DmNV_mau","PmNV","HgNV","DhNV","GpSGHV","MdSGHV","LbFV","AmFV","PoFV","LhFV","PcFV","EfFV","CcFV1","CcFV2"]

#Add the known viral dna ORFS 
cluster2_tab=cluster_tab.loc[~(cluster_tab['Names'].str.contains("_\\+_") | ~(~cluster_tab['Names'].str.contains("_\\-_")))]
cluster2_tab['Names2']= cluster2_tab['Names']
cluster2_tab['Names2'] = cluster2_tab['Names2'].str.replace('|'.join(List_viruses_name), '',regex=True)
record_dict = SeqIO.to_dict(SeqIO.parse( DNA_virus_ORFs, "fasta"))

for _, sub_cluster in cluster2_tab.groupby('Cluster'):
	filename=Output_directory+sub_cluster['Cluster'].iloc[0]+".dna"
	with open(filename, "a") as output:
		for index, row in sub_cluster.iterrows():
			print(">",row['Names'],sep="",file=output)
			print(record_dict[row['Names']].seq,file=output)


#Create core genome fasta clusers : 

Cluster_core_genes=[]
for filename in os.listdir("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs/"):
    if filename.endswith(".aa"):
      List_in_cluster=[]
      for record in SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs/"+filename, "fasta"):
        #print(record.id)
        #record=re.sub(".*_","",record.id)
        for viruses in  List_viruses_name_without_outgroupp:
          if viruses in record.id:
            if viruses not in List_in_cluster:
              List_in_cluster.append(viruses)
      if len(List_in_cluster) >= len(List_viruses_name_without_outgroupp):
        print(filename)
      if set( List_viruses_name_without_outgroupp).issubset(List_in_cluster) == True:
        print(filename)
        Cluster_core_genes.append(filename)
        print(List_in_cluster)


for core_clusters in Cluster_core_genes:
	with open("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs_core/"+core_clusters,"w") as output:
		record_list=[]
		for record in SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs/"+core_clusters, "fasta"):
			for viruses in  List_viruses_name:
				if viruses in record.id:
					if viruses not in record_list:
						record_list.append(viruses)
						print('>',viruses,sep="",file=output)
						print(record.seq,file=output)


print("\n")
print("All clusters created within the directory : ", Output_directory)

