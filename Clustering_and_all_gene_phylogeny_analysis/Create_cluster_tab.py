import os
import argparse
import networkx as nx
import pandas as pd 
import numpy as np

# Print out a message when the program is initiated. traitrise
print('----------------------------------------------------------------------------------------------\n')
print('                        Create Binary and cluster file                          .\n')
print('----------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Create Binary and cluster file ')
parser.add_argument("-c", "--cluster", help="The cluster output file")
parser.add_argument("-blast", "--blast", help="The all vs all blast file")
parser.add_argument("-binary", "--binary", help="The output binary file")
args = parser.parse_args()


#Usage example python3 /beegfs/home/bguinet/these_scripts_2/Create_cluster_tab.py --blast /beegfs/data/bguinet/LbFV_family_project/Clustering/All_versus_All_Predicted_and_known_ORFs.m8 --cluster /beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab --binary /beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs_binary.tab
blast_file=args.blast
output_cluster_file=args.cluster
output_binary_file=args.binary

#blast_file="/beegfs/data/bguinet/LbFV_family_project/Clustering/All_vs_All_All_known_and_new_viral_ORFs_result.m8"

df=pd.read_csv(blast_file,sep="\t",header=None)
df.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']

#df=df.loc[df['qcov'].ge(0.15) | df['tcov'].ge(0.15) ]
df=df.loc[df['bits'].gt(40)]
# Create the graph from the dataframe
g = nx.Graph()

g.add_edges_from(df[['query','target']].itertuples(index=False))
new = list(nx.connected_components(g))

mapped =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new) for node in component}

cluster = pd.DataFrame({'Cluster': mapped.values(), 'Names':mapped.keys()})

#Add alignment informations 
cluster_with_ali_info=cluster.copy()
df2=df[['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']]
df2.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
df3=df.append(df2)
df3=df3.sort_values(by=['evalue','bits'], ascending=[True,False])
df3=df3.loc[~(df3['query'] == df3['target'])]
df3 = df3.drop_duplicates(subset=['query'], keep='first')
cluster_with_ali_info=cluster_with_ali_info.merge(df3,right_on="query",left_on="Names")
#Get number of loci within clusters 
cluster_with_ali_info['Species_name'] =  cluster_with_ali_info['Names'].str.replace(".*_","")
cluster_with_ali_info['Nb_species'] = cluster_with_ali_info.groupby(['Cluster'])['Species_name'].transform('nunique')
cluster_with_ali_info_filtred = cluster_with_ali_info.loc[cluster_with_ali_info['Nb_species'].ge(4)]
print("Max evalue :", cluster_with_ali_info_filtred['evalue'].max())
print("Min pident :", cluster_with_ali_info_filtred['pident'].min())
print("Min bits :", cluster_with_ali_info_filtred['bits'].min())

#cluster.to_csv(output_cluster_file,sep=";",index=False)


#Create cluster files only for cluster with at least 4 species loci 
import re 
from Bio import SeqIO
#add length and Shorten the names :
cluster['Length']='NA'
cluster['Names2']='NA'
for loci in cluster['Names'].unique():
  if '_orf' in loci :
    new_loci=re.sub('_.*','',loci)
  else: 
    new_loci=re.sub('.*_','',loci)
  cluster.loc[cluster['Names']==loci,'Names2']=new_loci
  subdf=df.loc[df['query']==loci]
  if len(subdf)>=1:
    loci_length= subdf['qlen'].iloc[0]
  else:
    subdf=df.loc[df['target']==loci]
    loci_length= subdf['tlen'].iloc[0]
  cluster.loc[cluster['Names']==loci,'Length']=loci_length

#Count number of duplicates species loci within each clusters in order to find possibly wrongly associated clusters 
# set up group
g = cluster.groupby('Cluster')
# get unique values
cluster['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
cluster['Nb_duplicated'] = cluster['Number_unique_Names2'] - non_dup


cluster['Ratio_duplicated_unique'] = cluster['Nb_duplicated']/cluster['Number_unique_Names2']
duplicated_clusters=cluster.loc[cluster['Ratio_duplicated_unique'].gt(0.70) & cluster['Nb_duplicated'].ge(4)]

# Redifine more stringent clusters for those sequences 

df_bis=pd.read_csv(blast_file,sep="\t",header=None)
df_bis.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
df_bis=df_bis.loc[df_bis['qcov'].ge(0.25) | df_bis['tcov'].ge(0.25) ]
df_bis=df_bis.loc[df_bis['bits'].gt(47)]
# Create the graph from the dataframe
g_bis = nx.Graph()
g_bis.add_edges_from(df_bis[['query','target']].itertuples(index=False))
new_bis = list(nx.connected_components(g_bis))
mapped_bis =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new_bis) for node in component}
cluster_bis = pd.DataFrame({'Cluster': mapped_bis.values(), 'Names':mapped_bis.keys()})

#Add alignment informations 
cluster_with_ali_info_bis=cluster_bis.copy()
df2_bis=df_bis[['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']]
df2_bis.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
df3_bis=df_bis.append(df2_bis)
df3_bis=df3_bis.sort_values(by=['evalue','bits'], ascending=[True,False])
df3_bis=df3_bis.loc[~(df3_bis['query'] == df3_bis['target'])]
df3_bis = df3_bis.drop_duplicates(subset=['query'], keep='first')
cluster_with_ali_info_bis=cluster_with_ali_info_bis.merge(df3_bis,right_on="query",left_on="Names")
#Get number of loci within clusters 
cluster_with_ali_info_bis['Species_name'] =  cluster_with_ali_info_bis['Names'].str.replace(".*_","")
cluster_with_ali_info_bis['Nb_species'] = cluster_with_ali_info_bis.groupby(['Cluster'])['Species_name'].transform('nunique')
cluster_with_ali_info_filtred_bis = cluster_with_ali_info_bis.loc[cluster_with_ali_info_bis['Nb_species'].ge(4)]
print("Max evalue :", cluster_with_ali_info_filtred_bis['evalue'].max())
print("Min pident :", cluster_with_ali_info_filtred_bis['pident'].min())
print("Min bits :", cluster_with_ali_info_filtred_bis['bits'].min())

# Keep only cluster with too many duplicated sequences 
cluster_bis=cluster_bis.loc[cluster_bis['Names'].isin(duplicated_clusters['Names'].unique())]

#Remove from the original clusters the redifined clusters 
cluster=cluster.loc[~cluster['Names'].isin(cluster_bis['Names'].unique())]
cluster_bis['Cluster']=cluster_bis['Cluster']+'_redefined'
cluster=cluster.append(cluster_bis)

import re 
from Bio import SeqIO
#add length and Shorten the names :
cluster['Length']='NA'
cluster['Names2']='NA'
for loci in cluster['Names'].unique():
  if '_orf' in loci :
    new_loci=re.sub('_.*','',loci)
  else: 
    new_loci=re.sub('.*_','',loci)
  cluster.loc[cluster['Names']==loci,'Names2']=new_loci
  subdf=df.loc[df['query']==loci]
  if len(subdf)>=1:
    loci_length= subdf['qlen'].iloc[0]
  else:
    subdf=df.loc[df['target']==loci]
    loci_length= subdf['tlen'].iloc[0]
  cluster.loc[cluster['Names']==loci,'Length']=loci_length

#cluster.to_csv(output_cluster_file,sep=";",index=False)

#cluster=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab",sep=";")

#Print clusters containing onl LbFV-likes sequences 

list_LbFV_like=["LbFV","PoFV","LhFV","PcFV","EfFV","CcFV1","CcFV2"]

cluster[cluster.groupby('Cluster')['Names2'].transform(lambda x : set(x) == set(list_LbFV_like))]


# cluster = pd.read_csv(output_cluster_file,sep=";")

def to_dict_remove_dups(sequences):
      return {record.id: record for record in sequences}

record_dict = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Clustering/All_known_and_new_viral_ORFs.faa","fasta"))



#COrect pif-3 according to annie 
cluster.loc[cluster['Names'].str.contains("CcFV1_orf028"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79975_LbFV")]['Cluster'].iloc[0]

#Mannually add WSSV  into right clusters

cluster.loc[cluster['Names'].str.contains("YP_009220649"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79978")]['Cluster'].iloc[0] #DNA_pol
cluster.loc[cluster['Names'].str.contains("YP_009220510"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79983")]['Cluster'].iloc[0] #PIF0
cluster.loc[cluster['Names'].str.contains("YP_009220545"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79952")]['Cluster'].iloc[0] #PIF1
cluster.loc[cluster['Names'].str.contains("YP_009220486"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79972")]['Cluster'].iloc[0] #PIF2
cluster.loc[cluster['Names'].str.contains("YP_009220581"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79975")]['Cluster'].iloc[0] #PIF3
cluster.loc[cluster['Names'].str.contains("YP_009220481"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79963")]['Cluster'].iloc[0] #PIF5
cluster.loc[cluster['Names'].str.contains("YP_009220590"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ80017")]['Cluster'].iloc[0] #P33

#Manually add AmFV into right clusters

cluster.loc[cluster['Names'].str.contains("AKY03143"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79978")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03146"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79983")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03129"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79952")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03169"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79972")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03157"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79975")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03126"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79963")]['Cluster'].iloc[0]

cluster.to_csv(output_cluster_file,sep=";",index=False)

cluster=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab",sep=";")

#cluster=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab",sep=";")

#Write the cluster files 
grouped = cluster.groupby('Cluster')
for group_name, group in grouped:
	subcluster=group.sort_values(by='Length',ascending=False)
	subcluster=subcluster.drop_duplicates(subset='Names2', keep="first")
	if len(subcluster)>1:
		with open ("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs/"+group_name+".aa","w") as output: 
			for index, row in subcluster.iterrows():
				if 'CoBV' in row['Names2']:
					continue
				else:
					print('>',row['Names2'],sep="",file=output)
					print(re.sub("\\*","",str(record_dict[row['Names']].seq)),file=output)

#Write the Core cluster files
list_LbFV_like=["LbFV","PoFV","LhFV","PcFV","EfFV","CcFV1","CcFV2"]

grouped = cluster.groupby('Cluster')
for group_name, group in grouped:
	subcluster=group.sort_values(by='Length',ascending=False)
	subcluster=subcluster.drop_duplicates(subset='Names2', keep="first")
	n=0
	for i in list_LbFV_like:
		if  i in subcluster['Names2'].unique():
			n+=1
		else:
			continue
	if n == len(list_LbFV_like):
		with open ("/beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_seqs_core/"+group_name+".aa","w") as output: 
			for index, row in subcluster.iterrows():
				if 'CoBV' in row['Names2']:
                                        continue 
				else:
					print('>',row['Names2'],sep="",file=output)
					print(re.sub("\\*","",str(record_dict[row['Names']].seq)),file=output)

print("Cluster file written to : ",output_cluster_file)

#Create binary table 

List_viruses_name=["AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","HzNV-1","HzNV-2","GbNV","OrNV","ToNV","DiNV","DmNV_kal","DmNV_tom","DmNV_esp","DmNV_mau","PmNV","HgNV","DhNV","GpSGHV","MdSGHV","LbFV","AmFV","PoFV","LhFV","PcFV","EfFV","DFV","CcFV1","CcFV2"]

#cluster="ALL_Predicted_and_known_ORFs_cluster.tab"
#cluster=pd.read_csv(cluster,header=0,sep="\t")
#cluster['Names2'] = cluster['Names'].str.extract(f'({"|".join(List_viruses_name)})', expand=False)

#Load LbFV ORFs number table 

ORF_table = pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Clustering/LbFV_protein_ORF_names.txt",sep="\t")

cluster= cluster.merge(ORF_table,on='Names',how='outer')

cluster['Names2']=cluster['Names2'].replace('LbFV', np.nan)

cluster.Names2.fillna(cluster.New, inplace=True)

cluster_save=cluster.copy()


#Add ORF numbers 
#list_ORF_names=[]
#for viruses in  ["PoFV","LhFV","PcFV","EfFV","DFV","CcFV1","CcFV2"]:
#	subOrf_df=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+viruses+"_free/Prodigal_prediction/Prodigal_"+viruses+"_free_prediction", sep="\t",comment="#", header=None)
#	subOrf_df.columns=['Scaffold','prog','type','start','end','score','strand','zero','metadata']
#	print(viruses)
#	n=1
#	for index, row in subOrf_df.iterrows():
#		list_ORF_names.append({'ORF_name1':"ORF"+str(n)+":"+row['Scaffold']+"_"+str(row['start'])+"_"+str(row['end'])+"_"+row['strand']+"_"+viruses, 'ORF_name2':row['Scaffold']+"_"+str(row['start'])+"_"+str(row['end'])+"_"+row['strand']+"_"+viruses})
#		n+=1


#cluster_save=cluster_save.merge(pd.DataFrame(list_ORF_names),right_on="ORF_name2", left_on="Names",how="left")
#cluster_save = cluster_save.drop_duplicates(subset=['Cluster', 'ORF_name1'], keep='first')



del cluster['New']
del cluster['Names']

List_ORFs=list(ORF_table['New'])

m = cluster['Names2'].isin(List_ORFs)

f = lambda x: [y for y in x if y]

cluster['ORFs_values'] = cluster['Cluster'].map(cluster['Names2'].where(m, '').groupby(cluster['Cluster']).agg(f))

cluster = cluster.explode('ORFs_values')


cluster=(pd.get_dummies(cluster.dropna(subset=['ORFs_values']), 
                     columns=['Names2'], 
                     prefix='', 
                     prefix_sep='')
        .groupby(['ORFs_values','Cluster'], as_index=False)
        .max())

cluster = cluster.loc[:, ~cluster.columns.str.startswith('LbFV_')]

cluster['ORF_numbers'] = cluster['ORFs_values'].str.replace('.*_ORF','').astype(int)
cluster.sort_values(by='ORF_numbers', ascending=True,inplace=True)

cluster=cluster[['ORFs_values', 'Cluster', 'LhFV', 'DFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','CpV', 'CuniNPV', 'DhNV', 'tom', 'GbNV', 'GpSGHV', 'HzNV-1', 'LdMNPV','MdSGHV', 'NeseNPV', 'OrNV',  'PmNV',  'ToNV','AcMNPV','AmFV']]


#Add ORF names 


cluster_save['ORF_name1'] = "NA"

i=0
for index, row in cluster_save.iterrows():
	New_name="NA"
	try:
		New_name=str(cluster_save['Names'].str.split("_orf")[i][0]) + "_orf" +  str(re.sub("_.*","",str(cluster_save['Names'].str.split("_orf")[i][1])))
	except:
		New_name=str(cluster_save['Names'].str.split("_orf")[i][0]) 
	cluster_save['ORF_name1'][i]=New_name
	i+=1

#cluster_save['ORF_name1'] = cluster_save['ORF_name1'].str.strip()
cluster_save.ORF_name1.fillna(cluster_save.Names, inplace=True)
cluster_save = cluster_save.pivot_table(index='Cluster', columns='Names2', values='ORF_name1', aggfunc='|'.join, fill_value=0)
cluster[['ORFs_values','Cluster']].join(cluster_save, on='Cluster').fillna(0)

cluster = cluster[['ORFs_values','Cluster']].join(cluster_save, on='Cluster').fillna(0)

cluster=cluster[['ORFs_values', 'Cluster', 'LhFV', 'DFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','CpV', 'CuniNPV', 'DhNV', 'tom', 'GbNV', 'GpSGHV', 'HzNV-1', 'LdMNPV','MdSGHV', 'NeseNPV', 'OrNV',  'PmNV',  'ToNV','AcMNPV','AmFV']]


cluster.to_csv(output_binary_file,sep=";",index=False)

print("Binary file written to : ",output_binary_file)
#cluster[['ORFs_values', 'Cluster', 'LhFV', 'DFV','EfFV','PoFV','PcFV','CcFV1','CcFV2']]

