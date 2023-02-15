import pandas as pd 
import numpy as np
import argparse
import re,os
import subprocess
from taxadb.taxid import TaxID
from Bio import SeqIO 

Output_table=pd.DataFrame(columns=['Taxid', 'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'Assembly', 'new_start', 'new_end', 'strand', 'Taxid_x', 'clade', 'class', 'cohort', 'family', 'genus', 'infraclass', 'infraorder', 'kingdom', 'no rank', 'order', 'parvorder', 'phylum', 'species', 'subclass', 'subfamily', 'suborder', 'subphylum', 'superfamily', 'superkingdom', 'superorder', 'Taxid_y'])

for file in os.listdir("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/"):
	if file.endswith("_result_filtred.m8"):
		tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/"+file,sep=";")
		Output_table=Output_table.append(tab)

for file in os.listdir("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_BIPAA_LBBE/"):
	if file.endswith("_result_filtred.m8"):
		tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_BIPAA_LBBE/"+file,sep=";")
		Output_table=Output_table.append(tab)
#

Output_table.reset_index(drop=True, inplace=True)


Output_table=Output_table.loc[~Output_table['target'].eq("to_remove")]

Output_table=Output_table.drop_duplicates()

#Remove overlapping sequences 
#is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
#Output_table['Overlapp_group'] = Output_table.sort_values(['Assembly','target', 'tstart', 'tend']) \
#                .groupby(['Assembly','target'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()


#Output_table A= Output_table.sort_values(['Assembly', 'target', 'evalue'], ascending=[True, True, True]) \
#        .groupby(Output_table['Overlapp_group']).head(1)


Output_table['count_duplicate'] = Output_table.groupby(['Assembly','target','query'])['target'].transform('size')
Putative_HSPs=Output_table.loc[Output_table['count_duplicate'].ge(2)]
Putative_HSPs=Putative_HSPs.loc[~Putative_HSPs['Assembly'].str.contains("Encarsia")]
#Putative_HSPs = Putative_HSPs.sort_values(['Assembly', 'target', 'query','evalue'], ascending=[True, True, True,True]) \
#        .groupby(Putative_HSPs['Overlapp_group']).head(1)
Putative_HSPs=Putative_HSPs.sort_values(['Assembly', 'target','query'], ascending=[True, True,True])
Putative_HSPs['diff_length_target']=np.nan
Putative_HSPs['diff_length_query']=np.nan
Putative_HSPs['HSP_group']=np.nan
#Putative_HSPs=Putative_HSPs.loc[Putative_HSPs['target'].str.contains("HG993125.1")]



#is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
#Putative_HSPs['Overlapp_HSP_query_group'] = Putative_HSPs.sort_values(['query','target','HSP_group', 'qstart', 'qend']) \
#                .groupby(['query','target','HSP_group'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

#Output_table A= Putative_HSPs.sort_values(['query', 'target', 'Overlapp_HSP_query_group','bits'], ascending=[True, True, True,False]) \
#        .groupby(Putative_HSPs['Overlapp_HSP_query_group']).head(1)

#Putative_HSPs['Count_Overlapp_HSP_query_group']=Putative_HSPs.groupby("Overlapp_HSP_query_group")["Overlapp_HSP_query_group"].transform("count")

#grouped = Putative_HSPs.groupby(['Assembly', 'query','target'])
#Putative_HSPs = Putative_HSPs.sort_values(['Assembly', 'target', 'query','evalue'], ascending=[True, True, True,True]) \
#        .groupby(Putative_HSPs['Overlapp_group']).head(1)

#Putative_HSPs=Putative_HSPs.loc[Putative_HSPs['target'].str.contains("HG993125.1")]



Putative_HSPs['HSP_group']=Putative_HSPs.groupby(['Assembly', 'query','target']).ngroup()
Putative_HSPs['HSP_min_target']=Putative_HSPs.groupby('HSP_group')["tend"].transform("min")
Putative_HSPs['HSP_max_target']=Putative_HSPs.groupby('HSP_group')["tstart"].transform("max")
Putative_HSPs['HSP_min_query']=Putative_HSPs.groupby('HSP_group')["qend"].transform("min")
Putative_HSPs['HSP_max_query']=Putative_HSPs.groupby('HSP_group')["qstart"].transform("max")
Putative_HSPs['diff_length_target']=abs(Putative_HSPs['HSP_min_target']-Putative_HSPs['HSP_max_target'])
Putative_HSPs['diff_length_query']=abs(Putative_HSPs['HSP_min_query']-Putative_HSPs['HSP_max_query'])
"""
number=1
for group_name, group in grouped:
  diff_target=abs(min(group['tend'])-max(group['tstart']))
  diff_query=abs(min(group['qend'])-max(group['qstart']))
  #min(group['qend'])
  #max(group['qstart'])
  Group=number
  for index, row in group.iterrows():
    #Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) & (Putative_HSPs['Assembly']==row['Assembly'])][['query','target','tstart','tend','HSP_group','count_duplicate','Overlapp_group','diff_length','qstart','qend']]
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) & (Putative_HSPs['Assembly']==row['Assembly']),"diff_length_target"]= diff_target
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) & (Putative_HSPs['Assembly']==row['Assembly']),"diff_length_query"]= diff_query
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) & (Putative_HSPs['Assembly']==row['Assembly']),"HSP_group"]= Group
    #Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) & (Putative_HSPs['Assembly']==row['Assembly'])][['query','target','tstart','tend','HSP_group','count_duplicate','Overlapp_group','diff_length','qstart','qend']]
    number+=1
"""



is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
Putative_HSPs['Overlapp_HSP_target_group'] = Putative_HSPs.sort_values(['query','target','HSP_group', 'tstart', 'tend']) \
                .groupby(['query','target','HSP_group'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

#Merge overlapping targets 
#grouped = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])


Putative_HSPs["tstart"] = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])["tstart"].transform("min")
Putative_HSPs["tend"] = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])["tend"].transform("max")
Putative_HSPs["qstart"] = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])["qstart"].transform("min")
Putative_HSPs["qend"] = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])["qend"].transform("max")
Putative_HSPs["bits"] = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])['bits'].transform('mean') #ok
Putative_HSPs["qcov"] = Putative_HSPs.groupby(['Assembly','query', 'target','Overlapp_HSP_target_group'])['qcov'].transform('sum')
"""
for group_name, group in grouped:
  new_tstart=min(group['tstart'])
  #Gene_name=row['Gene_name']
  new_tend=max(group['tend'])
  new_bit=sum(group['bits'])/len(group)
  new_qcov=sum(group['qcov'])
  new_qstart=min(group['qstart'])
  new_qend=max(group['qend'])
  #tot_length=new_tend-new_tstart
  #diff_query=row['diff_length_query']
  #diff_target=group['diff_length_target']
  #Group=row['HSP_group']
  Overlapp_HSP_target_group=group['Overlapp_HSP_target_group'].iloc[0]
  for index, row in group.iterrows():
    Putative_HSPs.loc[(Putative_HSPs['Overlapp_HSP_target_group']==Overlapp_HSP_target_group) &   (Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"tstart"]= new_tstart
    Putative_HSPs.loc[(Putative_HSPs['Overlapp_HSP_target_group']==Overlapp_HSP_target_group) & (Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"tend"]= new_tend
    Putative_HSPs.loc[(Putative_HSPs['Overlapp_HSP_target_group']==Overlapp_HSP_target_group) & (Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"qcov"]=new_qcov
    Putative_HSPs.loc[(Putative_HSPs['Overlapp_HSP_target_group']==Overlapp_HSP_target_group) & (Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"bits"]=new_bit
    Putative_HSPs.loc[(Putative_HSPs['Overlapp_HSP_target_group']==Overlapp_HSP_target_group) & (Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"qstart"]=new_qstart
    Putative_HSPs.loc[(Putative_HSPs['Overlapp_HSP_target_group']==Overlapp_HSP_target_group) & (Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"qend"]=new_qend

"""
#Remove  duplicates target overlapping 
Putative_HSPs = Putative_HSPs.drop_duplicates(subset=['Assembly','target','query','Overlapp_HSP_target_group'], keep='first')


Putative_HSPs['HSP_min_target']=Putative_HSPs.groupby('HSP_group')["tend"].transform("min")
Putative_HSPs['HSP_max_target']=Putative_HSPs.groupby('HSP_group')["tstart"].transform("max")
Putative_HSPs['HSP_min_query']=Putative_HSPs.groupby('HSP_group')["qend"].transform("min")
Putative_HSPs['HSP_max_query']=Putative_HSPs.groupby('HSP_group')["qstart"].transform("max")
Putative_HSPs['diff_length_target']=abs(Putative_HSPs['HSP_max_target']-Putative_HSPs['HSP_min_target'])
Putative_HSPs['diff_length_query']=abs(Putative_HSPs['HSP_max_query']-Putative_HSPs['HSP_min_query'])

#Group close tstart-tend coordinates < 5OObp away from each other to merge them
Putative_HSPs=Putative_HSPs.sort_values(['Assembly','target','query','HSP_group', 'tstart'], ascending=[True, True,True,True,True])



def get_group(g):
    s = g['tend'].shift().rsub(g['tstart']).gt(500)
    return s.cumsum().add(1-s.iloc[0]).astype(str).radd('G')

Putative_HSPs['HSP_close_target_groups'] = (Putative_HSPs
   .groupby(['Assembly','target','query','HSP_group'], group_keys=False)
   .apply(get_group)
)

# Do not merge if the overlap of the query is to high > 50bp 

is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
Putative_HSPs['Overlapp_HSP_query_group'] = Putative_HSPs.sort_values(['query','target','HSP_group','HSP_close_target_groups', 'qstart', 'qend']) \
                .groupby(['query','target','HSP_group','HSP_close_target_groups'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()


#Merge close targets 

Putative_HSPs["tstart"] = Putative_HSPs.groupby(['Assembly','query', 'target','HSP_close_target_groups'])["tstart"].transform("min")
Putative_HSPs["tend"] = Putative_HSPs.groupby(['Assembly','query', 'target','HSP_close_target_groups'])["tend"].transform("max")
Putative_HSPs["qstart"] = Putative_HSPs.groupby(['Assembly','query', 'target','HSP_close_target_groups'])["qstart"].transform("min")
Putative_HSPs["qend"] = Putative_HSPs.groupby(['Assembly','query', 'target','HSP_close_target_groups'])["qend"].transform("max")
Putative_HSPs["bits"] = Putative_HSPs.groupby(['Assembly','query', 'target','HSP_close_target_groups'])['bits'].transform('mean') #ok
Putative_HSPs["qcov"] = Putative_HSPs.groupby(['Assembly','query', 'target','HSP_close_target_groups'])['qcov'].transform('sum')
Putative_HSPs = Putative_HSPs.drop_duplicates(subset=['Assembly','target','query','HSP_close_target_groups'], keep='first')

####
#Check if hit overlapp within the query between putative HSP that could then be duplicate genes 
is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
Putative_HSPs['Overlapp_HSP_query_group'] = Putative_HSPs.sort_values(['query','target','HSP_group', 'qstart', 'qend']) \
                .groupby(['query','target','HSP_group'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()
#
Putative_HSPs['count_Overlapp_HSP_query_group']=Putative_HSPs.groupby("Overlapp_HSP_query_group")["Overlapp_HSP_query_group"].transform("count")



#Remove all processed HSPs from the previous tab 
Output_table=Output_table.loc[Output_table['count_duplicate'].lt(2)]

Output_table2=Output_table.append(Putative_HSPs)
"""
#Putative_HSPs.loc[Putative_HSPs['count_HSP'].gt(2),"diff_length"]= 70

# If the distance between the different HSPs is > 80, do not merge them 

Putative_HSPs=Putative_HSPs.loc[Putative_HSPs['diff_length_target'].lt(10000) & Putative_HSPs['diff_length_query'].gt(-30) & Putative_HSPs['count_HSP'].eq(2) | Putative_HSPs['diff_length_target'].lt(20000) & Putative_HSPs['diff_length_query'].gt(-30) & Putative_HSPs['count_HSP'].gt(2)]
Putative_HSPs=Putative_HSPs.loc[~Putative_HSPs['Count_Overlapp_HSP_query_group'].ge(2)]
Output_table['diff_length']=np.nan
Output_table['HSP_group']=np.nan
Putative_HSPs['new_start']=np.nan
Putative_HSPs['new_end']=np.nan
Putative_HSPs['strand']="+"
#Putative_HSPs['Gene_name']=Putative_HSPs['query'].str.replace("^.*_([^_]*).*$", "\\1")
Putative_HSPs.loc[Putative_HSPs['tstart']> Putative_HSPs['tend'],"strand"]="-"
m = Putative_HSPs['strand'].eq('-')
Putative_HSPs[['tstart','tend']] = np.where(m.to_numpy()[:, None], 
                               Putative_HSPs[['tend','tstart']], 
                               Putative_HSPs[['tstart','tend']])
#

grouped = Putative_HSPs.groupby(['Assembly','query', 'target'])
tot_n=len(grouped)
n=0
#Output_table['keep']=np.nan
Output_table['Gene_name']=Output_table['query'].str.replace("^.*_([^_]*).*$", "\\1")
#Merge sequences together 
for group_name, group in grouped:
  new_tstart=min(group['tstart'])
  #Gene_name=row['Gene_name']
  new_tend=max(group['tend'])
  new_bit=sum(group['bits'])/len(group)
  new_qcov=sum(group['qcov'])
  tot_length=new_tend-new_tstart
  diff_query=row['diff_length_query']
  diff_target=row['diff_length_target']
  Group=row['HSP_group']
  #Output_table.loc[(Output_table['Gene_name']==Gene_name) &  (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"keep"]="no"
  for index, row in group.iterrows():
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"new_tstart"]= new_tstart
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"new_tend"]= new_tend
    #Output_table.loc[(Output_table['Gene_name']==Gene_name) &  (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"keep"]="no"
    #Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"keep"]= "yes"
    Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"diff_length_query"]= diff_query
    Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"diff_length_target"]= diff_target
    Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"HSP_group"]= Group
    Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"qcov"]= new_qcov
    Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"mean_HSP_bit"]=new_bit
    Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"Tot_length"]=tot_length
    if row['strand']=="+":
      Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"tstart"]= new_tstart
      Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"tend"]= new_tend
    if row['strand']=="-":
      Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"tend"]= new_tstart
      Output_table.loc[(Output_table['query']==row['query']) & (Output_table['target']==row['target']) &(Putative_HSPs['Assembly']==row['Assembly']),"tstart"]= new_tend
  n+=1
  print(n,"/",tot_n)






Output_table['mean_HSP_bit'] = np.where(Output_table['mean_HSP_bit'].isna(), Output_table['bits'], Output_table['mean_HSP_bit'])
Output_table['Tot_length'] = np.where(Output_table['Tot_length'].isna(), Output_table['tend'] - Output_table['tstart'] , Output_table['Tot_length'])
"""

Output_table2 = Output_table2.reset_index()
Output_table2['Gene_name']=Output_table2['query'].str.replace("^.*_([^_]*).*$", "\\1")
#del Output_table['level_0']
del Output_table2['index']

Output_table2.loc[Output_table2['Assembly'].str.contains("GCA_018492975.1") & Output_table2['Gene_name'].str.contains("pif3"),"query"]="CcFV1_pif3"

List_filamentous_loci=[]
for loci in Output_table2['query']:
           for filamentous in ['LbFV','LhFV','DFV','PcFV','PoFV','CcFV1','CcFV2','EfFV']:
             if filamentous in loci:
               List_filamentous_loci.append(loci)

m = Output_table2['query'].isin(List_filamentous_loci)


Output_table2.loc[Output_table2['Assembly'].str.contains("GCA_018492975.1") & Output_table2['Gene_name'].str.contains("pif3"),"query"]="CuniNPV_pif3"


Output_table2 = Output_table2.reset_index()

Output_table2['evalue']=Output_table2['evalue'].replace([0.0,], 2.747000e-310)

Output_table2['ln_evalue'] = Output_table2['evalue'].apply(np.log10)

Output_table2['diff_evalue'] = (Output_table2.groupby(['Assembly','target','Gene_name'])['ln_evalue']
                       .transform(lambda d: d.mask(m).min() -  d.where(m).min() ))


#Keep the very best hit 
Output_table2['Overlapp_group'] = Output_table2.sort_values(['Assembly','target','Gene_name','tstart','tend']) \
                .groupby(['Assembly','target','Gene_name'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Output_table2['Tot_length'] = abs(Output_table2['tend'] - Output_table2['tstart'])
Output_table2 = Output_table2.sort_values(['Assembly', 'target', 'Tot_length'], ascending=[True, True, False]) \
        .groupby(Output_table2['Overlapp_group']).head(1)

#Define thresholds based on known cases 
# Select the mean diff_evalue among controls to be our threshold
Output_table2.loc[Output_table2['Assembly'].eq("GCA_001855655.1"),'species']="Leptopilina clavipes"

List_controls=["Leptopilina heterotoma","Leptopilina boulardi","Leptopilina clavipes"]

mean_control_diff_evalue=np.mean(Output_table2.loc[(Output_table2['species'].isin(List_controls)) & abs(Output_table2['diff_evalue'].gt(0)) & ~(Output_table2['query'].str.contains('dnapol|ATPase')) & ~(Output_table2['diff_evalue'].isna() | Output_table2['diff_evalue'].gt(1000000))]['diff_evalue'].unique())
min_control_diff_evalue=np.min(Output_table2.loc[(Output_table2['species'].isin(List_controls)) & Output_table2['diff_evalue'].gt(0) & ~(Output_table2['query'].str.contains('dnapol|ATPase')) & ~(Output_table2['diff_evalue'].isna() | Output_table2['diff_evalue'].gt(1000000))]['diff_evalue'].unique())


Output_table2.loc[Output_table2['diff_evalue'].isna(),"diff_evalue"]=100000000000
Output_table2=Output_table2.loc[Output_table2['diff_evalue'].ge(min_control_diff_evalue)]

# Output_table2[['Assembly','query','target','evalue','ln_evalue','diff_evalue','species']]


# keep only if best hit is filamentous 
Output_table2.sort_values(by=["evalue"], inplace = True)
Output_table2.drop_duplicates(subset =['Assembly','target','Gene_name'],keep = "first", inplace = True) 
Output_table2.loc[Output_table2['Assembly'].str.contains("GCA_018492975.1") & Output_table2['Gene_name'].str.contains("pif3"),"query"]="CcFV1_pif3"
List_filamentous_loci=[]
for loci in Output_table2['query']:
           for filamentous in ['LbFV','LhFV','DFV','PcFV','PoFV','CcFV1','CcFV2','EfFV']:
             if filamentous in loci:
               List_filamentous_loci.append(loci)

Output_table2=Output_table2.loc[Output_table2['query'].isin(List_filamentous_loci)]
Output_table2.loc[Output_table2['Assembly'].str.contains("GCA_018492975.1") & Output_table2['Gene_name'].str.contains("pif3"),"query"]="CuniNPV_pif3"
 # Filter 
print("Summary table")`
Output_table2=Output_table2.loc[~ (Output_table2['query'].str.contains("AmFV"))]
Output_table2=Output_table2.loc[~ (Output_table2['query'].str.contains("ATPase|helicase") &  Output_table2['evalue'].gt(0.0000000005))]
g = Output_table2.groupby(['Assembly'])
Output_table2['Nb_EVEs'] = g['query'].transform('nunique') 


# Add sum of tlen to filter small cumulative scaffold that could correspond to Free-living viruses 
Output_table2['Cumulative_scaffold_length'] = Output_table2.groupby(['Assembly'])['tlen'].transform('sum')


#Manually add missing taxonomy data 

Missing_taxo_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Missing_taxonomy.txt",sep="\t")

for assembly in Missing_taxo_tab['Assembly'].unique():
        Output_table2.loc[Output_table2['Assembly'].eq(assembly),'species'] = Missing_taxo_tab.loc[Missing_taxo_tab['Assembly'].eq(assembly)]['species'].iloc[0]
        Output_table2.loc[Output_table2['Assembly'].eq(assembly),'family'] = Missing_taxo_tab.loc[Missing_taxo_tab['Assembly'].eq(assembly)]['family'].iloc[0]
        Output_table2.loc[Output_table2['Assembly'].eq(assembly),'order'] = Missing_taxo_tab.loc[Missing_taxo_tab['Assembly'].eq(assembly)]['order'].iloc[0]


# Do another step to add informations about genomes 


New_tax_tab=pd.DataFrame(columns=['Taxid','Assembly', 'Taxid_x', 'clade', 'class', 'cohort', 'family', 'genus', 'infraclass', 'infraorder', 'kingdom', 'no rank', 'order', 'parvorder', 'phylum', 'species', 'subclass', 'subfamily', 'suborder', 'subphylum', 'superfamily', 'superkingdom', 'superorder', 'Taxid_y'])
#New_tax_tab['Assembly']=Output_table2.loc[Output_table2['species'].isna()]['Assembly'].unique()


n_assembly=len(Output_table2.loc[Output_table2['species'].isna()]['Assembly'].unique())
n=0
for assembly_name in Output_table2.loc[Output_table2['species'].isna()]['Assembly'].unique():
	taxid=subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/esearch -db assembly -q '"+re.sub("_result","",assembly_name)+"' | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid", shell=True,capture_output=True, text=True)
	taxid=taxid.stdout
	taxid=re.sub(".*\t",'',taxid)
	taxid=re.sub("\n.*",'',taxid)
	taxids = TaxID(dbtype='sqlite', dbname='/beegfs/data/bguinet/taxadb2/taxadb.sqlite')
	Taxid_dictionnary = {}
	Taxid_dictionnary[taxid]= taxids.lineage_name(taxid,ranks=True,reverse=True)
	New_tax_tab.loc[New_tax_tab['Assembly'].eq(assembly_name),'Taxid']=taxid
	try:
		taxid_db=pd.DataFrame({k: dict(v) for k, v in Taxid_dictionnary.items()}).T
		taxid_db['Taxid']=taxid
		taxid_db['Assembly']=assembly_name
		#New_tax_tab=New_tax_tab.merge(taxid_db, left_on='Taxid', right_index=True,how='outer')
		New_tax_tab=New_tax_tab.append(taxid_db)
	except:
		print("did not work")
	print(n, " / ", n_assembly)
	n+=1
        print(New_tax_tab)


subOutput_table2= Output_table2.loc[Output_table2['species'].isna()]
subOutput_table2=subOutput_table2[['index', 'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'Assembly', 'new_start', 'new_end', 'strand', 'Overlapp_group', 'count_duplicate', 'HSP_group', 'Gene_name', 'ln_evalue', 'diff_evalue', 'Nb_EVEs', 'Cumulative_scaffold_length']]
subOutput_table2=subOutput_table2.merge(New_tax_tab,on="Assembly",how='outer')

subOutput_table3=Output_table2.loc[~Output_table2['species'].isna()]

Output_table2=subOutput_table2.append(subOutput_table3)
Output_table2['species'] = np.where(df['Y']>=50, 'yes', 'no')

#Manually add missing taxonomy data 

#Missing_taxo_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Missing_taxonomy.txt",sep="\t")

#for assembly in Missing_taxo_tab['Assembly'].unique():#
#	Output_table2.loc[Output_table2['Assembly'].eq(assembly),'species'] = Missing_taxo_tab.loc[Missing_taxo_tab['Assembly'].eq(assembly)]['species'].iloc[0]
#        Output_table2.loc[Output_table2['Assembly'].eq(assembly),'family'] = Missing_taxo_tab.loc[Missing_taxo_tab['Assembly'].eq(assembly)]['family'].iloc[0]
#        Output_table2.loc[Output_table2['Assembly'].eq(assembly),'order'] = Missing_taxo_tab.loc[Missing_taxo_tab['Assembly'].eq(assembly)]['order'].iloc[0]

#Output_table2.loc[Output_table2['Assembly']=="GCA_946863875.1","species"] = "Alloplasta piceator"
#Output_table2.loc[Output_table2['Assembly']=="GCA_946863875.1","family"] = "Ichneumonidea"
#Output_table2.loc[Output_table2['Assembly']=="GCA_946863875.1","order"] = "Hymenoptera"

#Output_table2.loc[Output_table2['Assembly']=="GCA_020882685.1","species"] = "Aphelinus atriplicis"
#Output_table2.loc[Output_table2['Assembly']=="GCA_020882685.1","family"] = "Aphelinidae"
#Output_table2.loc[Output_table2['Assembly']=="GCA_020882685.1","order"] = "Hymenoptera"

#Output_table2.loc[Output_table2['Assembly']=="GCA_018237165.1","species"] = "Hypophylla argenissa"
#Output_table2.loc[Output_table2['Assembly']=="GCA_018237165.1","family"] = "Riodinidae"
#Output_table2.loc[Output_table2['Assembly']=="GCA_018237165.1","order"]= "Lepidoptera"

#Output_table2.loc[Output_table2['Assembly']=="GCA_022816925.1","species"] = "Trichacis sp"
#Output_table2.loc[Output_table2['Assembly']=="GCA_022816925.1","family"] = "Platygastridae"
#Output_table2.loc[Output_table2['Assembly']=="GCA_022816925.1","order"] = "Hymenoptera"

Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";",index=False)

Output_table2=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")

#Output_table2.loc[Output_table2['target'].eq("JAAXZA010261437.1") & Output_table2['query'].str.contains("PcFV_LbFVorf99"),"strand"] = "-"
#Output_table2.loc[Output_table2['target'].eq("JAAXZA010154453.1") & Output_table2['query'].str.contains("DFV_pif3"),"strand"] = "-"
#Output_table2.loc[Output_table2['target'].eq("JAAXZA010217671.1") & Output_table2['query'].str.contains("CcFV2_LbFVorf99"),"tend"]=4741

Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";",index=False)

Output_table2=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")
# Manually remove loci not from Filamentous thanks to the phylogeny checked by hand
#nonfilamentous_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Table_with_nonfilamentous_loci.txt",sep="\t")

#Output_table2['full_name']=Output_table2['species']+":"+Output_table2['target']+":"+Output_table2['tstart'].astype(str)
#Output_table2.loc[Output_table2['target'].str.contains("DWUI010000029.1"),"full_name"]="nan:DWUI010000029.1:3546"
#Output_table2.loc[Output_table2['target'].str.contains("DWUI010000269"),"full_name"]="nan:DWUI010000269.1:803"


#Output_table2=Output_table2.loc[~Output_table2['full_name'].isin(nonfilamentous_table['full_name'])]

# Generate binary table for the heatmap plot 
#Binary_Output_table2=Output_table2.copy()
#Binary_Output_table2 = Binary_Output_table2.drop_duplicates(subset = ["Assembly", "target","Gene_name"])
#Binary_Output_table2 = Binary_Output_table2.join(pd.crosstab(Output_table2['Assembly'], Output_table2['Gene_name']).add_prefix(''), on='Assembly')
#Binary_Output_table2=Binary_Output_table2[['Assembly','genus','species','family','order','Cumulative_scaffold_length','Nb_EVEs','38K', 'ATPase', 'Integrase2', 'LbFVorf102', 'LbFVorf20', 'LbFVorf23', 'LbFVorf5', 'LbFVorf87', 'LbFVorf92', 'LbFVorf99', 'PDDEXK', 'ac38', 'ac81', 'dnapol', 'helicase', 'lcat', 'lef4', 'lef8', 'lef9', 'p33', 'p74', 'pif1', 'pif2', 'pif3', 'pif5']]
#Binary_Output_table2=Binary_Output_table2.drop_duplicates()
#Binary_Output_table2.sort_values(by=["order","family","genus"], inplace = True)
#save 
#Binary_Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes/Binary_filamentous_heatmap_table_EVEs_FL.txt",sep=";",index=False)

#Download assemblies if they do not exist 
n=0
number=len(Output_table2["Assembly"].unique())
for assembly in Output_table2["Assembly"].unique():
	if os.path.exists(assembly+".fna"):
		n+=1
		#print("############")
		#print(n ,"/", number)
	else:
		print(assembly)
		print("find . -name '"+assembly+"*.fna*' -type f -delete")
		print("esearch -db assembly -query "+assembly+" </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | while read -r url ; do fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ; wget $url/$fname ; done")
		print("rename 's@"+assembly+".*@"+assembly+".fna.gz@gs' *"+assembly+"*.fna.gz")
		print("gzip -d "+assembly +".fna.gz")
		subprocess.run("find . -name '"+assembly+"*.fna*' -type f -delete",shell=True)
		subprocess.run("esearch -db assembly -query "+assembly+" </dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | while read -r url ; do fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ; wget $url/$fname ; done", shell=True)
		subprocess.run("rename 's@"+assembly+".*@"+assembly+".fna.gz@gs' *"+assembly+"*.fna.gz",shell=True)
		subprocess.run("gzip -d "+assembly +".fna.gz",shell=True)
		n+=1
		print("############")
		print(n ,"/", number)

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

All_naldaviricetes=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/All_filamentoviridae_core_genes.aa", "fasta"))

# Generate cluster seq files 

#First create file and add Naldaviricetes genes 

#Output_table2=Output_table3
#Output_table3=Output_table2.copy()


#Output_table2=Output_table2.loc[Output_table2['query'].str.contains("orf94")]
for gene in Output_table2['Gene_name'].unique():
                print("###", gene)
                with open("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/"+gene+".aa","w") as output:
                        sub_Output_table2=Output_table2.loc[Output_table2['Gene_name'].eq(gene)]
                        #Add naldaviricetes 
                        for key, value  in All_naldaviricetes.items():
                                if gene in key:
                                        print(">",All_naldaviricetes[key].id,sep="",file=output)
                                        print(All_naldaviricetes[key].seq,file=output)
# Fore each euk species, add the EVEs 

Added_loci=[]
Added_assembly=[]
for assembly in Output_table2['Assembly'].unique():
	print("####", assembly)
	if assembly in Added_assembly:
		print(assembly ,' already added')
	else:
		try:
			assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/"+assembly+".fna", "fasta"))
		except:
			assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_BIPAA_LBBE/"+assembly+".fna", "fasta"))
		for gene in Output_table2.loc[Output_table2['Assembly'].eq(assembly)]['Gene_name'].unique():
			print("###", gene)
			with open("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/"+gene+".aa","a") as output:
				sub_Output_table2=Output_table2.loc[Output_table2['Gene_name'].eq(gene) & Output_table2['Assembly'].eq(assembly)]
				for index, row in sub_Output_table2.iterrows():
					if row['strand'] =="+":
							loci_name=re.sub(" ","_",row['species'])+"-"+row['target']+":"+str(row['tstart'])+"-"+str(row['tend'])+"("+row['strand']+")|"+row['order']
							if loci_name in Added_loci:
								print("loci already added")
							else:
								print(">",loci_name,sep="",file=output)
								print(str(assembly_genome[row['target']].seq[row['tstart']-1:row['tend']].translate()),file=output)
								Added_loci.append(loci_name)
							#print(">",row['species'],"-",row['target'],":",int(row['tstart']),"-",int(row['tend']),"(",row['strand'],")|",row['order'],sep="")
							#print(str(assembly_genome[row['target']].seq[int(row['tstart'])-1:int(row['tend'])].translate()))
					if row['strand'] =="-":
							loci_name=re.sub(" ","_",row['species'])+"-"+row['target']+":"+str(row['tstart'])+"-"+str(row['tend'])+"("+row['strand']+")|"+row['order']
							if loci_name in Added_loci:
								print("loci already added")
							else:
								print(">",loci_name,sep="",file=output)
								print(str(assembly_genome[row['target']].seq[row['tstart']-1:row['tend']].reverse_complement().translate()),file=output)
								Added_loci.append(loci_name)
						#print(">",row['species'],"-",row['target'],":",int(row['tstart']),"-",int(row['tend']),"(",row['strand'],")|",row['order'],sep="")
						#print(str(assembly_genome[row['target']].seq[int(row['tstart'])-1:int(row['tend'])].reverse_complement().translate()))
	Added_assembly.append(assembly)



#Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";",index=False)
Output_table2=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")


# Manually remove loci not from Filamentous thanks to the phylogeny checked by hand
#nonfilamentous_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Table_with_nonfilamentous_loci.txt",sep="\t")

#Output_table2['full_name']=Output_table2['species']+":"+Output_table2['target']+":"+Output_table2['tstart'].astype(int).astype(str)
#Output_table2.loc[Output_table2['target'].str.contains("DWUI010000029.1"),"full_name"]="nan:DWUI010000029.1:3546"
#Output_table2.loc[Output_table2['target'].str.contains("DWUI010000269"),"full_name"]="nan:DWUI010000269.1:803"

#Output_table2=Output_table2.loc[~Output_table2['full_name'].isin(nonfilamentous_table['full_name'])]

#####################################################
# Find ORFs along the scaffolds with candidate loci #
#####################################################

# Write the scaffolds containing the hits into a file for each assembly 
tot=len(Output_table2['Assembly'].unique())
n=0

Loaded_scaffolds=[]
if os.path.exists("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.dna") :
	Loaded_scaffolds_dic= to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.dna","fasta"))
	Loaded_scaffolds=list(Loaded_scaffolds_dic.keys())

with open ("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.dna","a") as output:
	for assembly in Output_table2.loc[~Output_table2['Assembly_target'].isin(Loaded_scaffolds)]['Assembly'].unique():
		subOutput_table2=Output_table2.loc[Output_table2['Assembly'].eq(assembly)]
		subOutput_table2['Assembly_target']=subOutput_table2['Assembly']+":"+subOutput_table2['target']
		subOutput_table2=subOutput_table2.loc[~subOutput_table2['Assembly_target'].isin(Loaded_scaffolds)]
		try:
			assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/"+assembly+".fna", "fasta"))
		except:
			assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_BIPAA_LBBE/"+assembly+".fna", "fasta"))
		for scaffold in subOutput_table2['target'].unique():
			#print(">",assembly,":",assembly_genome[scaffold].id,sep="",file=output)
			assembly_genome[scaffold].id=assembly+":"+assembly_genome[scaffold].id
			assembly_genome[scaffold].description=assembly+":"+assembly_genome[scaffold].description
			#print(textwrap.fill(str(assembly_genome[scaffold].seq),width=80),file=output)
			SeqIO.write(assembly_genome[scaffold],output,"fasta")
		n+=1
		print(n,"/",tot)


All_scaffolds=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.dna", "fasta"))


with open("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds2.dna", "w") as handle:
    SeqIO.write(All_scaffolds.values(), handle, "fasta") 
# Write the hits 
subprocess.run("cat /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/*.aa >> /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa" , shell=True)

# Run python ORF finder 
import pathlib  

Candidate_scaffolds="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds2.dna"
Candidate_scaffolds_ORFs_bed="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.bed"
Candidate_scaffolds_ORFs_dna="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.fna"
Candidate_scaffolds_ORFs_aa="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.aa"
outdir="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/"

subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/orfipy " +Candidate_scaffolds+"  --ignore-case  --procs 20 --outdir "+outdir+" --bed "+Candidate_scaffolds_ORFs_bed+" --min 150 --start ATG --dna "+Candidate_scaffolds_ORFs_dna +" --pep "+Candidate_scaffolds_ORFs_aa, shell=True)

subprocess.run(" sed -i 's@).*@)@g' " +Candidate_scaffolds_ORFs_aa, shell=True)
subprocess.run("sed -i 's@ \[@;@g' " +Candidate_scaffolds_ORFs_aa, shell=True)
subprocess.run("sed -i 's@\]@@g' " +Candidate_scaffolds_ORFs_aa, shell=True)
# Extract and translate the ORFs
ORF_bed=pd.read_csv(Candidate_scaffolds_ORFs_bed,sep="\t",header=None)
ORF_bed.columns=['Scaffold_name','ORF_start','ORF_end','ORF_name','zero','ORF_strand']

ORF_bed['ORF_name2']=ORF_bed['ORF_name'].str.replace(";.*","")
ORF_bed['ORF_name2']=ORF_bed['ORF_name2'].str.replace("ID=","")
ORF_bed['ORF_name']=ORF_bed['ORF_name'].str.replace("ID=","")
ORF_bed['dict_name']=ORF_bed['ORF_name2']+";"+ORF_bed['ORF_start'].astype(str)+"-"+ORF_bed['ORF_end'].astype(str)+"("+ORF_bed['ORF_strand']+")"
#Only keep ORFs around 5000bp from candidats to note keep to much ORFs to blast 
Selected_ORF_name=[]
n_tot=len(Output_table2['Assembly_target'].unique())
n=0

for assembly_target in Output_table2['Assembly_target'].unique():
	subORF_bed=ORF_bed.loc[ORF_bed['Scaffold_name'].eq(assembly_target)]
	subOutput_table2=Output_table2.loc[Output_table2['Assembly_target'].eq(assembly_target)]
	for index, row in subOutput_table2.iterrows():
		min_start=row['tstart']-5000
		max_start=row['tend']+5000
		for loci in subORF_bed.loc[subORF_bed['ORF_start'].lt(max_start) &  subORF_bed['ORF_end'].gt(min_start)]['dict_name'].unique():
			Selected_ORF_name.append(loci)
	n+=1
	print(n, "/", n_tot)

Candidate_scaffolds_ORFs_aa_filtred="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds_filtrer.aa"

ORF_dict=SeqIO.to_dict(SeqIO.parse(Candidate_scaffolds_ORFs_aa, "fasta"))

ORF_dict2=[ORF_dict.pop(x) for x in Selected_ORF_name if x in ORF_dict.keys()]

with open(Candidate_scaffolds_ORFs_aa_filtred,'w') as output:
	for loci in Selected_ORF_name:
		print(">",ORF_dict[loci].id,sep="",file=output)
		print(ORF_dict[loci].seq,file=output)

"""
Loaded_ORF=[]
with open(Candidate_scaffolds_ORFs_aa,"w") as output:
        record_dict=to_dict_remove_dups(SeqIO.parse(Candidate_scaffolds_ORFs_dna,"fasta"))
        for species in ORF_bed['Scaffold_name'].unique():
                for index, row in ORF_bed.loc[ORF_bed['Scaffold_name'].str.contains(species)].iterrows():
			ORF_name=row['ORF_name2']+';'+str(row['ORF_start'])+"-"+str(row['ORF_end'])+"("+row['ORF_strand']+")"
			if ORF_name in Loaded_ORF:
				continue
			else:
                        	print(">",ORF_name,sep="",file=output)
	                        print(record_dict[row['ORF_name2']].seq.translate(),file=output)
				Loaded_ORF.append(ORF_name)
"""

# Run mmseqs between ORFs and previously selected candidates

Filtred_loci="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa"
subprocess.run("sed -i 's@ @_@g'  /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa", shell=True)

Filtred_loci_db=outdir+"All_candidate_hits_db"

Filtred_loci_vs_Predicted_ORFs_result=outdir+"Filtred_loci_vs_Predicted_ORFs_result"
Filtred_loci_vs_Predicted_ORFs_temp=outdir+"Filtred_loci_vs_Predicted_ORFs_tpm"
Filtred_loci_vs_Predicted_ORFs_table=outdir+"Filtred_loci_vs_Predicted_ORFs_result.m8"

#Filtred_loci_vs_Predicted_ORFs_result1=outdir+"Filtred_loci_vs_Predicted_ORFs_result1"
#Filtred_loci_vs_Predicted_ORFs_temp1=outdir+"Filtred_loci_vs_Predicted_ORFs_tpm1"
#Filtred_loci_vs_Predicted_ORFs_table1=outdir+"Filtred_loci_vs_Predicted_ORFs_result1.m8"

#Filtred_loci_vs_Predicted_ORFs_result2=outdir+"Filtred_loci_vs_Predicted_ORFs_result2"
#Filtred_loci_vs_Predicted_ORFs_temp2=outdir+"Filtred_loci_vs_Predicted_ORFs_tpm2"
#Filtred_loci_vs_Predicted_ORFs_table2=outdir+"Filtred_loci_vs_Predicted_ORFs_result2.m8"



#Split file into 3 files 


#subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/pyfasta split -n 3 "+ Candidate_scaffolds_ORFs_aa, shell=True)


Candidate_scaffolds_ORFs_aa="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.aa"
#Candidate_scaffolds_ORFs_aa1="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.aa.1"
#Candidate_scaffolds_ORFs_aa2="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.aa.2"

Predicted_ORFs_db=outdir+"scaffold_orfipy_db"
#Predicted_ORFs_db1=outdir+"scaffold_orfipy_db1"
#Predicted_ORFs_db2=outdir+"scaffold_orfipy_db2"

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Candidate_scaffolds_ORFs_aa_filtred + " "+ Predicted_ORFs_db , shell=True)
#subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Candidate_scaffolds_ORFs_aa1 + " "+ Predicted_ORFs_db1 , shell=True)
#subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Candidate_scaffolds_ORFs_aa2 + " "+ Predicted_ORFs_db2 , shell=True)

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Filtred_loci + " "+ Filtred_loci_db , shell=True)

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  Filtred_loci_db  +" "+ Predicted_ORFs_db + " "+  Filtred_loci_vs_Predicted_ORFs_result + " "+ Filtred_loci_vs_Predicted_ORFs_temp + " -e 0.000001 --threads 30", shell=True)
#subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  Filtred_loci_db  +" "+ Predicted_ORFs_db1 + " "+  Filtred_loci_vs_Predicted_ORFs_result1 + " "+ Filtred_loci_vs_Predicted_ORFs_temp1 + " -e 0.000001 --threads 30", shell=True)
#subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  Filtred_loci_db  +" "+ Predicted_ORFs_db2 + " "+  Filtred_loci_vs_Predicted_ORFs_result2 + " "+ Filtred_loci_vs_Predicted_ORFs_temp2 + " -e 0.000001 --threads 30", shell=True)


subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_loci_db +" "+ Predicted_ORFs_db + " "+  Filtred_loci_vs_Predicted_ORFs_result + " "+ Filtred_loci_vs_Predicted_ORFs_table + " --threads 10", shell=True)
#subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_loci_db +" "+ Predicted_ORFs_db1 + " "+  Filtred_loci_vs_Predicted_ORFs_result1 + " "+ Filtred_loci_vs_Predicted_ORFs_table1 + " --threads 10", shell=True)
#subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_loci_db +" "+ Predicted_ORFs_db2 + " "+  Filtred_loci_vs_Predicted_ORFs_result2 + " "+ Filtred_loci_vs_Predicted_ORFs_table2 + " --threads 10", shell=True)

Filtred_loci_vs_Predicted_ORFs_table_ALL=outdir+"Filtred_loci_vs_Predicted_ORFs_result.m8"

#Merge all files 

subprocess.run(" cat "+ Filtred_loci_vs_Predicted_ORFs_table0 + " " + Filtred_loci_vs_Predicted_ORFs_table1 + " " + Filtred_loci_vs_Predicted_ORFs_table0 + ">> "+ Filtred_loci_vs_Predicted_ORFs_table_ALL , shell=True)

ORF_vs_EVE_table=pd.read_csv(Filtred_loci_vs_Predicted_ORFs_table_ALL,sep="\t",header=None)
ORF_vs_EVE_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov','qcov']
# Keep only self matching hits 
ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['bits'].ge(50)]


#Output_table2 = pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")

ORF_vs_EVE_table['query_EVE_species']=ORF_vs_EVE_table['query'].str.replace(":.*","")
#ORF_vs_EVE_table['query_EVE_species'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV"),ORF_vs_EVE_table['query'].str.replace(".*_",""),inplace=True)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(\\+\\)","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(-\\)","")
ORF_vs_EVE_table['ORF_strand']="NA"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("\\+\\)"),'ORF_strand'] ="+"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("-\\)"),'ORF_strand'] ="-"
#ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_end']=ORF_vs_EVE_table['ORF_start'].str.replace(".*-","").astype(int)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("-.*","").astype(int)
ORF_vs_EVE_table['target_ORF_species']=ORF_vs_EVE_table['target'].str.replace("_ORF.*","")
ORF_vs_EVE_table['target_ORF_species1']=ORF_vs_EVE_table['target_ORF_species'].str.replace(".*:","")
ORF_vs_EVE_table['Assembly']=ORF_vs_EVE_table['target_ORF_species'].str.replace(":.*","")
ORF_vs_EVE_table=ORF_vs_EVE_table.merge(Output_table2[['Assembly','species']], on ="Assembly",how="left")
ORF_vs_EVE_table['target_ORF_species']=ORF_vs_EVE_table['species'].str.replace(" ","_")+"-"+ORF_vs_EVE_table['target_ORF_species1']

ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query'].str.replace(":.*","")
ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query_EVE_scaffold'].str.replace(".*-","")
#ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query'].str.replace("-.*",""),inplace=True)
#ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query_EVE_scaffold'].str.split('_').str[:-1].str.join('_'),inplace=True)
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target'].str.replace(".*:","")
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target_ORF_scaffold'].str.replace("_ORF.*","")

ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query_EVE_species'] == ORF_vs_EVE_table['target_ORF_species']]

# Remove duplicate 
ORF_vs_EVE_table = ORF_vs_EVE_table.drop_duplicates()

# find overlapping ORFs within the same scaffold
ORF_vs_EVE_table.rename(columns={'query': 'ORF_query',
                   'target': 'ORF_target'},
          inplace=True, errors='raise')

Output_table2['Names']=Output_table2['species'].str.replace(" ","_")+"-"+Output_table2['target']+":"+Output_table2['tstart'].astype(int).astype(str)+"-"+Output_table2['tend'].astype(int).astype(str)+"("+Output_table2['strand']+")|"+Output_table2['order']
Output_table2_ORFs=Output_table2.merge(ORF_vs_EVE_table[['ORF_query','ORF_target','ORF_start','ORF_end','ORF_strand']],left_on="Names",right_on="ORF_query",how="left")


Output_table2_ORFs['Overlapp_ORF_EVEs']= np.nan
Output_table2_ORFs['ORF_perc']= np.nan



for index, row in Output_table2_ORFs.loc[Output_table2_ORFs['ORF_start'].ge(0)].iterrows():
        overlapp="no"
        if (int(row['tstart']) <= int(row['ORF_start'])) & (int(row['tend']) <= int(row['ORF_end'])) & (int(row['tend']) > int(row['ORF_start'])):
                #print("left")
                #print(row['tstart']," : ", row['tend'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-left"
        elif (int(row['tstart']) >= int(row['ORF_start'])) &  (int(row['tend']) >= int(row['ORF_end'])) & (int(row['tstart']) < int(row['ORF_end'])):
                #print("right")
                #print(row['tstart']," : ", row['tend'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-right"
        elif (int(row['tstart']) >= int(row['ORF_start'])) &  (int(row['tend']) <= int(row['ORF_end'])):
                overlapp = "yes-inside"
        elif (int(row['tstart']) <= int(row['ORF_start'])) & (int(row['tend']) >= int(row['ORF_end'])) :
                #print("outside")
                #print(row['tstart']," : ", row['tend'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-outside"
        else:
                overlapp = "no"
        Output_table2_ORFs.loc[Output_table2_ORFs['ORF_query'].eq(row['ORF_query']) & Output_table2_ORFs['ORF_end'].eq(row['ORF_end']) & Output_table2_ORFs['tend'].eq(row['tend']),"Overlapp_ORF_EVEs"]=overlapp
        if overlapp =="no":
                continue
        else:
                Output_table2_ORFs.loc[Output_table2_ORFs['ORF_query'].eq(row['ORF_query']) & Output_table2_ORFs['ORF_end'].eq(row['ORF_end']) & Output_table2_ORFs['tend'].eq(row['tend']),"ORF_perc"]= (int(row['ORF_end'])- int(row['ORF_start']))/ (int(row['tend'])- int(row['tstart']))
        #print("\n")

# remove non-overlapping ORFS with EVEs 
Output_table2_ORFs.loc[Output_table2_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_start"] = np.nan
Output_table2_ORFs.loc[Output_table2_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_end"] = np.nan

Output_table2_ORFs.loc[Output_table2_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & Output_table2_ORFs['ORF_perc'].lt(0.5),"ORF_start"] = np.nan
Output_table2_ORFs.loc[Output_table2_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & Output_table2_ORFs['ORF_perc'].lt(0.5),"ORF_end"] = np.nan

# find overlapping ORFs within the same scaffold

Output_table2_ORFs['full_name']=Output_table2_ORFs['species']+":"+Output_table2_ORFs['target']+":"+Output_table2_ORFs['tstart'].astype(int).astype(str)

is_overlapped = lambda x: x['ORF_start'] >= x['ORF_end'].shift(fill_value=-1)
Output_table2_ORFs['ORF_overlapp_group'] = Output_table2_ORFs.sort_values(['full_name', 'ORF_start', 'ORF_end']) \
                .groupby(['full_name'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Output_table2_ORFs.loc[Output_table2_ORFs['ORF_start'].isna(),"ORF_overlapp_group"] = np.nan

# Calculate ORF_perc_best_hit 
Output_table2_ORFs['Best_hit_ORF_perc']= ((Output_table2_ORFs['ORF_end']-Output_table2_ORFs['ORF_start'])/3)/Output_table2_ORFs['qlen']

Output_table2_ORFs.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results_ORFs.tab",sep=";",index=False)


Output_table2_ORFs=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results_ORFs.tab",sep=";")


# Manually remove loci not from Filamentous thanks to the phylogeny checked by hand
nonfilamentous_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Table_with_nonfilamentous_loci.txt",sep="\t")

Output_table2_ORFs['full_name']=Output_table2_ORFs['species'].str.replace(" ","_")+"-"+Output_table2_ORFs['target']+"_"+Output_table2_ORFs['tstart'].astype(str)+"-"+Output_table2_ORFs['tend'].astype(str)+"("+Output_table2_ORFs['strand']+")|"+Output_table2_ORFs['order']

#To check correspondance
number_removed_loci=0
for loci in nonfilamentous_table['full_name']:
	if len(Output_table2_ORFs.loc[Output_table2_ORFs['full_name'].eq(loci)])==0:
		print(loci)
	else:
		number_removed_loci+=1

print("Number removed loci : ",number_removed_loci) 
#Number removed loci :  373

Output_table2_ORFs=Output_table2_ORFs.loc[~Output_table2_ORFs['full_name'].isin(nonfilamentous_table['full_name'].unique())]


#Output_table2_ORFs=Output_table2_ORFs.loc[Output_table2_ORFs['diff_evalue'].gt(min_control_diff_evalue)]

Output_table2_ORFs_sum=Output_table2_ORFs.drop_duplicates(subset = ["Assembly","target"])
Output_table2_ORFs_sum['Cumulative_scaffold_length'] = Output_table2_ORFs_sum.groupby(['Assembly'])['tlen'].transform('sum')
del Output_table2_ORFs['Cumulative_scaffold_length']
Output_table2_ORFs=Output_table2_ORFs.merge(Output_table2_ORFs_sum[['Assembly','Cumulative_scaffold_length']],on="Assembly")

Output_table2_ORFs = Output_table2_ORFs.drop_duplicates()

# Add pseudogenized informations 
Output_table2_ORFs['pseudogenized']="no"

Output_table2_ORFs_grouped = Output_table2_ORFs.groupby('Gene_name')
#Output_table2_ORFs['Assembly_target']= Output_table2_ORFs['Assembly']+":"+ Output_table2_ORFs['target']
Output_table2_ORFs['full_name']=Output_table2_ORFs['species'].str.replace(" ","_")+"-"+Output_table2_ORFs['target']+":"+Output_table2_ORFs['tstart'].astype(str)+"-"+Output_table2_ORFs['tend'].astype(str)+"("+Output_table2_ORFs['strand']+")|"+Output_table2_ORFs['order']
# iterate over each group
for group_name, df_group in Output_table2_ORFs_grouped:
	subprocess.run("sed -i 's@ @_@g'  /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/"+group_name+".aa", shell=True)
	fasta_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/"+group_name+".aa", "fasta"))
	for index, row in df_group.iterrows():
		loci_name=re.sub(" ","_",row['species'])+"-"+row['target']+":"+str(row['tstart'])+"-"+str(row['tend'])+"("+row['strand']+")|"+row['order']
		if "*" in fasta_dict[loci_name].seq:
			Output_table2_ORFs.loc[Output_table2_ORFs['full_name'].eq(loci_name),"pseudogenized"]="yes"



#Remove some duplicate 

Output_table2_ORFs=Output_table2_ORFs.sort_values(['Assembly', 'ORF_start'], ascending=[True, False])
Output_table2_ORFs = Output_table2_ORFs.drop_duplicates(subset=['Assembly','target','Gene_name','tstart','tend'], keep='first')

Output_table2_ORFs['ORF_completness']=np.nan
Output_table2_ORFs['Assembly_gene_name']=Output_table2_ORFs['Assembly']+'-'+Output_table2_ORFs['Gene_name']

m = Output_table2_ORFs.Best_hit_ORF_perc > 0.75

s1 = Output_table2_ORFs.loc[m, 'Assembly_gene_name'].unique()
s2 = Output_table2_ORFs.loc[~m, 'Assembly_gene_name'].unique()

m1 = Output_table2_ORFs['Assembly_gene_name'].isin(s1)
m2 = Output_table2_ORFs['Assembly_gene_name'].isin(s2)

Output_table2_ORFs['ORF_completness'] = np.select([m1 & ~m2, ~m1 & m2], [0, 2], default=1)

Output_table2_ORFs=Output_table2_ORFs.sort_values(['Assembly', 'Gene_name'], ascending=[True, False])


both_completness_tab=Output_table2_ORFs.loc[(Output_table2_ORFs['ORF_completness'].eq(1))][['Assembly','species','target','query','evalue','tstart','tend','ORF_start','ORF_end','Best_hit_ORF_perc','ORF_completness','Gene_name','pseudogenized']]
both_completness_tab.drop_duplicates(subset=['Assembly','Gene_name'], keep='first')

Output_table2_ORFs.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results_ORFs_pseudogenized.tab",sep=";",index=False)
Output_table2_ORFs=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results_ORFs_pseudogenized.tab",sep=";")


#Add status 
Output_table2_ORFs['Status']=0
Output_table2_ORFs.loc[Output_table2_ORFs['Best_hit_ORF_perc'].ge(0.70) & Output_table2_ORFs['pseudogenized'].eq("no"),'Status']=1
Output_table2_ORFs.loc[Output_table2_ORFs['Best_hit_ORF_perc'].ge(0.70) & Output_table2_ORFs['pseudogenized'].eq("yes"),'Status']=2
Output_table2_ORFs.loc[Output_table2_ORFs['Best_hit_ORF_perc'].lt(0.70) & Output_table2_ORFs['pseudogenized'].eq("yes"),'Status']=3
Output_table2_ORFs.loc[Output_table2_ORFs['Best_hit_ORF_perc'].lt(0.70) & Output_table2_ORFs['pseudogenized'].eq("no"),'Status']=4
Output_table2_ORFs.loc[Output_table2_ORFs['Best_hit_ORF_perc'].isna() & Output_table2_ORFs['pseudogenized'].eq("yes"),'Status']=3
Output_table2_ORFs.loc[Output_table2_ORFs['Best_hit_ORF_perc'].isna() & Output_table2_ORFs['pseudogenized'].eq("no"),'Status']=4

# Generate binary table for the heatmap plot 

g = Output_table2_ORFs.groupby(['species', 'Gene_name'])
binary_table=g['Status'].agg(lambda x: '-'.join(x.astype('str').sort_values().unique())).unstack()
binary_table.index.name = 'species'
binary_table.reset_index(inplace=True)
binary_table=binary_table.fillna(0)

#Add order 

#binary_table=binary_table.merge(Output_table2_ORFs[['species','order']],on="species")


#Add higher cumulative length 
#binary_table=binary_table.merge(Output_table2_ORFs[['species','Cumulative_scaffold_length']].sort_values(by=["Cumulative_scaffold_length"],ascending=False).drop_duplicates(subset = "species"),on="species")
binary_table = binary_table.drop_duplicates()
binary_table = binary_table.merge(binary_table.set_index('species').apply(pd.Series.value_counts, axis=1).fillna(0).add_prefix('count_'), left_on='species', right_index=True)
binary_table['Nb_EVEs']=29-binary_table['count_0']
#Add order 
binary_table=binary_table.merge(Output_table2_ORFs[['species','order']],on="species")


#Add higher cumulative length 
binary_table=binary_table.merge(Output_table2_ORFs[['species','Cumulative_scaffold_length']].sort_values(by=["Cumulative_scaffold_length"],ascending=False).drop_duplicates(subset = "species"),on="species")



binary_table.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Binary_filamentous_heatmap_table_EVEs_FL3.txt",sep=";",index=False)
binary_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Binary_filamentous_heatmap_table_EVEs_FL2.txt",sep=";")

"""
#manualyy fix AC81 for clavipes 
Output_table2_ORFs.loc[Output_table2_ORFs['Gene_name'].str.contains("ac81") & Output_table2_ORFs['species'].str.contains("clavipes"),"Best_hit_ORF_perc"]=0.90
Binary_Output_table2=Output_table2_ORFs.copy()
Binary_Output_table2 = Binary_Output_table2.drop_duplicates(subset = ["Assembly", "target","Gene_name"])
Binary_Output_table2 = Binary_Output_table2.join(pd.crosstab(Output_table2_ORFs['Assembly'], Output_table2_ORFs['Gene_name']).add_prefix(''), on='Assembly')
Binary_Output_table2 = Binary_Output_table2.drop_duplicates()
Binary_Output_table2=Binary_Output_table2[['Assembly','full_name','tstart','tend','genus','species','family','order','Best_hit_ORF_perc','Cumulative_scaffold_length','Nb_EVEs','38K', 'ATPase', 'Integrase2', 'LbFVorf102', 'LbFVorf20', 'LbFVorf23', 'LbFVorf5', 'LbFVorf54','LbFVorf87', 'LbFVorf92','LbFVorf94', 'LbFVorf99', 'PDDEXK', 'ac38', 'ac81', 'dnapol', 'helicase', 'lcat', 'lef4', 'lef5','lef8', 'lef9', 'p33', 'p74', 'pif1', 'pif2', 'pif3', 'pif5']]
Binary_Output_table2=Binary_Output_table2.drop_duplicates(subset = "Assembly")
Binary_Output_table2.sort_values(by=["order","family","genus"], inplace = True)

# Add 2 if the ORF is not long enough compared to the best hit to be potentialy functional 

cols=['38K', 'ATPase', 'Integrase2', 'LbFVorf102', 'LbFVorf20', 'LbFVorf23', 'LbFVorf5','LbFVorf54', 'LbFVorf87', 'LbFVorf92','LbFVorf94', 'LbFVorf99', 'PDDEXK', 'ac38', 'ac81', 'dnapol', 'helicase', 'lcat', 'lef4','lef5', 'lef8', 'lef9', 'p33', 'p74', 'pif1', 'pif2', 'pif3', 'pif5']

for species in Binary_Output_table2['species'].unique():
	subOutput_table2_ORFs=Output_table2_ORFs.loc[Output_table2_ORFs['species'].eq(species)]
	for gene in cols:
		if len(subOutput_table2_ORFs.loc[subOutput_table2_ORFs['Gene_name'].eq(gene)])>0:
			if subOutput_table2_ORFs.loc[subOutput_table2_ORFs['Gene_name'].eq(gene)]['Best_hit_ORF_perc'].max() >= 0.70:
				for assembly in subOutput_table2_ORFs['Assembly'].unique():
					Binary_Output_table2.loc[Binary_Output_table2['Assembly'].eq(assembly),gene]=1
			else:
				for assembly in subOutput_table2_ORFs['Assembly'].unique():
					Binary_Output_table2.loc[Binary_Output_table2['Assembly'].eq(assembly),gene]=2


#print(subBinary_Output_table2[gene].iloc[0]
#mask = Binary_Output_table2[cols].ne(0) & Binary_Output_table2['Best_hit_ORF_perc'].gt(0.85).to_numpy()[:, None]
#Binary_Output_table2[cols] = np.where(mask, 1, 2) * Binary_Output_table2[cols].ne(0)

#save 
Binary_Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Binary_filamentous_heatmap_table_EVEs_FL.txt",sep=";",index=False)



#Output_table2_ORFs['Pseudogenized']="no"
#Filtred_loci_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa","fasta"))










# Table output 
#Output_table2=Output_table2.loc[ ~(Output_table2['Cumulative_scaffold_length'].lt(200000) & Output_table2['Nb_EVEs'].gt(13))]

print(Output_table2[['Assembly','species','family','order','Nb_EVEs','Cumulative_scaffold_length']].drop_duplicates(subset = "Assembly"))
print(len(Output_table2[['Assembly','species','family','order','Nb_EVEs']].drop_duplicates(subset = "Assembly")))

"""
