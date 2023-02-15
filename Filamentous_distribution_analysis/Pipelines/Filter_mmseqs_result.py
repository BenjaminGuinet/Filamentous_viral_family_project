import pandas as pd 
import numpy as np
import argparse
import re,os
import subprocess
from taxadb.taxid import TaxID
from Bio import SeqIO 

# Print out a message when the program is initiated.
print('----------------------------------------------------------------\n')
print('                        Filter output from Blast .\n')
print('----------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add strand information and change coordinate position in the blast output file')
parser.add_argument("-f", "--file", help=" the mmseqs blast file in .m8 format")
parser.add_argument("-o", "--output", help=" the desired output filtred mmseqs blast file")

args = parser.parse_args()

file =args.file
output= args.output

# Example usage 
#python3 Filter_mmseqs_result.py -f /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/GCA_905404225.1_result.m8  -o /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/GCA_905404225.1_result_filtred.m8

"""
if os.path.exists(output) :
 print ("File already exists, openin this file...")
 Output_table=pd.read_csv(output,sep=";")
else:
 print ("We have to create the All file... ")
 Output_table=pd.DataFrame(columns=['Assembly','query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov','new_start','new_end','strand','Taxid'])
"""

filename = re.sub(".*\\/","",file)

if filename.endswith("result.m8"):
    #print(filename)
      input=filename 
      assembly_name=re.sub(".m8","",input)
      #output=re.sub(".m8","_filtred.m8",input)
      print ("Parsing : ", filename)
      if os.stat(file).st_size != 0:
        Mmseqs_tab=pd.read_csv(file,sep="\t",header=None)
        #if  filename in ['GCA_019393585.1_result_orf94.m8','GCA_015476485.1_result_orf94.m8','GCA_011634795.1_result_orf94.m8']:
        #   Mmseqs_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','tstart','tend','qstart','qend','evalue','bits','qlen','tlen','tcov','qcov'] 
        Mmseqs_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
        print(filename)
        print(Mmseqs_tab)
        #Filter 
        Mmseqs_tab=Mmseqs_tab.loc[Mmseqs_tab['qcov'].ge(0.25)]
        Mmseqs_tab=Mmseqs_tab.loc[Mmseqs_tab['bits'].ge(50)]
        Mmseqs_tab['Assembly']=re.sub("_result","",assembly_name)
        if len(Mmseqs_tab)!=0:
         Mmseqs_tab['new_start']=np.nan
         Mmseqs_tab['new_end']=np.nan
         Mmseqs_tab['strand']="+"
         Mmseqs_tab.loc[Mmseqs_tab['tstart']> Mmseqs_tab['tend'],"strand"]="-"
         m = Mmseqs_tab['strand'].eq('-')
         Mmseqs_tab[['tstart','tend']] = np.where(m.to_numpy()[:, None], 
                                        Mmseqs_tab[['tend','tstart']], 
                                        Mmseqs_tab[['tstart','tend']])
         #
         #Group overlapping loci 
         is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
         #manual stuff 
         piece=Mmseqs_tab.iloc[0]
         piece['target']="to_remove"
         Mmseqs_tab=Mmseqs_tab.append(piece)
         Mmseqs_tab.reset_index(drop=True, inplace=True)
         #
         List_filamentous_loci=[]
         for loci in Mmseqs_tab['query']:
           for filamentous in ['LbFV','LhFV','DFV','PcFV','PoFV','CcFV1','CcFV2']:
             if filamentous in loci:
               List_filamentous_loci.append(loci)
         #
         sub_Mmseqs_tab=Mmseqs_tab.loc[Mmseqs_tab['query'].isin(List_filamentous_loci)]
         Mmseqs_tab=Mmseqs_tab.loc[Mmseqs_tab['Assembly'].isin(sub_Mmseqs_tab['Assembly'])]
         #
         #Load species informations about the assembly 
         if len(Mmseqs_tab[~Mmseqs_tab['query'].str.contains('ATPase')])!=0:
          #print("Number of filtred ORFs : ", len(Mmseqs_tab['query'].unique()))
          #print("Number of filtred loci : ", len(Mmseqs_tab))
          taxid=subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/esearch -db assembly -q '"+re.sub("_result","",assembly_name)+"' | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid", shell=True,capture_output=True, text=True)
          taxid=taxid.stdout
          taxid=re.sub(".*\t",'',taxid)
          taxid=re.sub("\n.*",'',taxid)
          taxids = TaxID(dbtype='sqlite', dbname='/beegfs/data/bguinet/taxadb2/taxadb.sqlite')
          Taxid_dictionnary = {}
          Taxid_dictionnary[taxid]= taxids.lineage_name(taxid,ranks=True,reverse=True)
          Mmseqs_tab['Taxid']=taxid
          try:
           taxid_db=pd.DataFrame({k: dict(v) for k, v in Taxid_dictionnary.items()}).T
           taxid_db['Taxid']=taxid
           Mmseqs_tab=Mmseqs_tab.merge(taxid_db, left_on='Taxid', right_index=True,how='outer')
          except:
           print("")
          Mmseqs_tab.to_csv(output,sep=";",index=False)


#if os.path.exists(output):
output="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes/ALL_results.m8"

Output_table3=Output_table.copy()

Output_table=pd.read_csv(output,sep=";")


Output_table=Output_table.append(Output_table3)
Output_table = Output_table.drop_duplicates(subset = ["Assembly", "target","query","tstart"])
#Output_table = Output_table.drop_duplicates()
# Add orf94 analysis 
output2="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes/ALL_results_orf94.m8"
Output_table2=pd.read_csv(output2,sep=";")

Output_table=Output_table.append(Output_table2)
Output_table.reset_index(drop=True, inplace=True)


Output_table=Output_table.loc[~Output_table['target'].eq("to_remove")]

#Remove overlapping sequences 
is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
Output_table['Overlapp_group'] = Output_table.sort_values(['Assembly','target', 'tstart', 'tend']) \
                .groupby(['Assembly','target'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

#Output_table = Output_table.sort_values(['Assembly', 'target', 'evalue'], ascending=[True, True, True]) \
#        .groupby(Output_table['Overlapp_group']).head(1)


Output_table = Output_table.reset_index()
Output_table['Gene_name']=Output_table['query'].str.replace("^.*_([^_]*).*$", "\\1")
#del Output_table['level_0']
del Output_table['index']

List_filamentous_loci=[]
for loci in Output_table['query']:
           for filamentous in ['LbFV','LhFV','DFV','PcFV','PoFV','CcFV1','CcFV2']:
             if filamentous in loci:
               List_filamentous_loci.append(loci)

m = Output_table['query'].isin(List_filamentous_loci)

Output_table = Output_table.reset_index()

Output_table['evalue']=Output_table['evalue'].replace([0.0,], 2.747000e-310)

Output_table['ln_evalue'] = Output_table['evalue'].apply(np.log10)

Output_table['diff_evalue'] = (Output_table.groupby(['Assembly','target','Gene_name'])['ln_evalue']
                       .transform(lambda d: d.mask(m).min() -  d.where(m).min() ))

#Define thresholds based on known cases 
# Select the mean diff_evalue among controls to be our threshold
List_controls=["Leptopilina heterotoma","Leptopilina boulardi","Leptopilina clavipes"]
mean_control_diff_evalue=np.mean(Output_table.loc[(Output_table['species'].isin(List_controls)) & Output_table['diff_evalue'].gt(0) & ~(Output_table['query'].str.contains('dnapol|ATPase')) & ~(Output_table['diff_evalue'].isna() | Output_table['diff_evalue'].gt(1000000))]['diff_evalue'].unique())
min_control_diff_evalue=np.min(Output_table.loc[(Output_table['species'].isin(List_controls)) & Output_table['diff_evalue'].gt(0) & ~(Output_table['query'].str.contains('dnapol|ATPase')) & ~(Output_table['diff_evalue'].isna() | Output_table['diff_evalue'].gt(1000000))]['diff_evalue'].unique())


Output_table.loc[Output_table['diff_evalue'].isna(),"diff_evalue"]=100000000000
Output_table2=Output_table.loc[Output_table['diff_evalue'].ge(min_control_diff_evalue)]

# Output_table2[['Assembly','query','target','evalue','ln_evalue','diff_evalue','species']]


# keep only if best hit is filamentous 
Output_table2.sort_values(by=["evalue"], inplace = True)
Output_table2.drop_duplicates(subset =['Assembly','target','Gene_name'],keep = "first", inplace = True) 
List_filamentous_loci=[]
for loci in Output_table2['query']:
           for filamentous in ['LbFV','LhFV','DFV','PcFV','PoFV','CcFV1','CcFV2']:
             if filamentous in loci:
               List_filamentous_loci.append(loci)

Output_table2=Output_table2.loc[Output_table2['query'].isin(List_filamentous_loci)]

 # Filter 

print("Summary table")`
Output_table2=Output_table2.loc[~ (Output_table2['query'].str.contains("AmFV"))]
Output_table2=Output_table2.loc[~ (Output_table2['query'].str.contains("ATPase|helicase") &  Output_table2['evalue'].gt(0.0000000005))]
g = Output_table2.groupby(['Assembly'])
Output_table2['Nb_EVEs'] = g['query'].transform('nunique') 

# Add sum of tlen to filter small cumulative scaffold that could correspond to Free-living viruses 
Output_table2['Cumulative_scaffold_length'] = Output_table2.groupby(['Assembly'])['tlen'].transform('sum')

#Manually add missing taxonomy data 

Output_table2.loc[Output_table2['Assembly']=="GCA_946863875.1","species"] = "Alloplasta piceator"
Output_table2.loc[Output_table2['Assembly']=="GCA_946863875.1","family"] = "Ichneumonidea"
Output_table2.loc[Output_table2['Assembly']=="GCA_946863875.1","order"] = "Hymenoptera"

Output_table2.loc[Output_table2['Assembly']=="GCA_020882685.1","species"] = "Aphelinus atriplicis"
Output_table2.loc[Output_table2['Assembly']=="GCA_020882685.1","family"] = "Aphelinidae"
Output_table2.loc[Output_table2['Assembly']=="GCA_020882685.1","order"] = "Hymenoptera"

Output_table2.loc[Output_table2['Assembly']=="GCA_018237165.1","species"] = "Hypophylla argenissa"
Output_table2.loc[Output_table2['Assembly']=="GCA_018237165.1","family"] = "Riodinidae"
Output_table2.loc[Output_table2['Assembly']=="GCA_018237165.1","order"]= "Lepidoptera"

Output_table2.loc[Output_table2['Assembly']=="GCA_022816925.1","species"] = "Trichacis sp"
Output_table2.loc[Output_table2['Assembly']=="GCA_022816925.1","family"] = "Platygastridae"
Output_table2.loc[Output_table2['Assembly']=="GCA_022816925.1","order"] = "Hymenoptera"

Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";",index=False)

Output_table2=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")

Output_table2.loc[Output_table2['target'].eq("JAAXZA010261437.1") & Output_table2['query'].str.contains("PcFV_LbFVorf99"),"strand"] = "-"
Output_table2.loc[Output_table2['target'].eq("JAAXZA010154453.1") & Output_table2['query'].str.contains("DFV_pif3"),"strand"] = "-"
Output_table2.loc[Output_table2['target'].eq("JAAXZA010217671.1") & Output_table2['query'].str.contains("CcFV2_LbFVorf99"),"tend"]=4741

Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";",index=False)

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




#Load NR results 
Output_table2=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")


NR_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/All_clusters_NR_result.m8",sep="\t",header=None)
NR_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','taxid','taxname','taxlineage']
NR_table=NR_table.loc[NR_table['bits'].gt(50)]

import numpy as np
NR_table['order']=np.nan
NR_table.loc[NR_table['taxlineage'].str.contains("acter"),"order"]="Bacteria"
NR_table.loc[NR_table['taxlineage'].str.contains("ukaryot"),"order"]="Eukaryota"
NR_table.loc[NR_table['taxlineage'].str.contains("irus"),"order"]="Virus"
NR_table.loc[NR_table['taxlineage'].str.contains("rchaea"),"order"]="Archaea"

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


Output_table2=Output_table2.loc[Output_table2['query'].str.contains("lcat")]
for gene in Output_table2['Gene_name'].unique():
                print("###", gene)
                with open("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters2/"+gene+".aa","w") as output:
                        sub_Output_table2=Output_table2.loc[Output_table2['Gene_name'].eq(gene)]
                        #Add naldaviricetes 
                        for key, value  in All_naldaviricetes.items():
                                if gene in key:
                                        print(">",All_naldaviricetes[key].id,sep="",file=output)
                                        print(All_naldaviricetes[key].seq,file=output)
list_non_viral=[]
list_non_viral2=[]
# Fore each euk species, add the EVEs 
for assembly in Output_table2['Assembly'].unique():
	print("####", assembly)
	#assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/"+assembly+".fna", "fasta"))
	for gene in Output_table2.loc[Output_table2['Assembly'].eq(assembly)]['Gene_name'].unique():
		print("###", gene)
		with open("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters2/"+gene+".aa","a") as output:
			sub_Output_table2=Output_table2.loc[Output_table2['Gene_name'].eq(gene) & Output_table2['Assembly'].eq(assembly)]
			for index, row in sub_Output_table2.iterrows():
				try:
					row_name=re.sub(" ","_",row['species'])+"-"+row['target']+":"+str(row['tstart'])+"-"+str(row['tend'])+"("+row['strand']+")|"+row['order']
					print(row_name)
				except:
					continue
				subNR_table=NR_table.loc[NR_table['query'].eq(row_name)]
				subNR_table = subNR_table.drop_duplicates(subset = "taxlineage")
				#Take the 3 best hits
				for i in list(subNR_table.loc[subNR_table['query'].eq(row_name) & (~subNR_table['order'].eq("Virus"))].iloc[0:3]['target']):
					if i in list_non_viral:
						continue
					else:
						if subNR_table.loc[subNR_table['target'].eq(i)]['taxname'].iloc[0] +":"+gene in list_non_viral2:
							continue
						else:
							list_non_viral2.append(subNR_table.loc[subNR_table['target'].eq(i)]['taxname'].iloc[0] +":"+gene)
							list_non_viral.append(i)
				#assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes/"+row['Assembly']+".fna", "fasta"))
				#assembly=row['Assembly']
				#row['start']=sub_Output_tabl.loc[sub_Output_table2['Assembly'].eq(assembly)]['tstart'].iloc[0]
				#end=sub_Output_table2.loc[sub_Output_table2['Assembly'].eq(assembly)]['tend'].iloc[0]
				#strand=sub_Output_table2.loc[sub_Output_table2['Assembly'].eq(assembly)]['strand'].iloc[0]
				#scaffold=sub_Output_table2.loc[sub_Output_table2['Assembly'].eq(assembly)]['target'].iloc[0]
				#species=re.sub(" ","_",sub_Output_table2.loc[sub_Output_table2['Assembly'].eq(assembly)]['species'].iloc[0])
				#order=sub_Output_table2.loc[sub_Output_table2['Assembly'].eq(assembly)]['order'].iloc[0]
				if row['strand'] =="+":
						print(">",row['species'],"-",row['target'],":",row['tstart'],"-",row['tend'],"(",row['strand'],")|",row['order'],sep="",file=output)
						print(str(assembly_genome[row['target']].seq[row['tstart']-1:row['tend']].translate()),file=output)
						#print(">",row['species'],"-",row['target'],":",int(row['tstart']),"-",int(row['tend']),"(",row['strand'],")|",row['order'],sep="")
						#print(str(assembly_genome[row['target']].seq[int(row['tstart'])-1:int(row['tend'])].translate()))
				if row['strand'] =="-":
						print(">",row['species'],"-",row['target'],":",row['tstart'],"-",row['tend'],"(",row['strand'],")|",row['order'],sep="",file=output)
						print(str(assembly_genome[row['target']].seq[row['tstart']-1:row['tend']].reverse_complement().translate()),file=output)
						#print(">",row['species'],"-",row['target'],":",int(row['tstart']),"-",int(row['tend']),"(",row['strand'],")|",row['order'],sep="")
						#print(str(assembly_genome[row['target']].seq[int(row['tstart'])-1:int(row['tend'])].reverse_complement().translate()))



#For each gene, add non-viral sequences
from Bio import Entrez
from Bio import SeqIO 
Entrez.email = "benjamin.guinet95@gmail.com"
Entrez.api_key = "30bf99cff0e43d6827934fa6ab127f3b5f09"
list_non_viral = list(dict.fromkeys(list_non_viral))


for gene in Output_table2['Gene_name'].unique():
	with open("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters2/"+gene+".aa","a") as output:
		subOutput_table2=Output_table2.loc[Output_table2['Gene_name'].eq(gene)]	
		subOutput_table2['query2']=subOutput_table2['species'].str.replace(" ","_")+"-"+subOutput_table2['target']+":"+subOutput_table2['tstart'].astype(str)+"-"+subOutput_table2['tend'].astype(str)+"("+subOutput_table2['strand']+")|"+subOutput_table2['order']
		#subOutput_table2['query2']=re.sub(" ","_",subOutput_table2['species'])+"-"+subOutput_table2['target']+":"+str(subOutput_table2['tstart'])+"-"+str(subOutput_table2['tend'])+"("+subOutput_table2['strand']+")|"+subOutput_table2['order']
		subNR_table=NR_table.loc[NR_table['query'].isin(subOutput_table2['query2'])]
		subNR_table=subNR_table.loc[subNR_table['target'].isin(list_non_viral)]
		for i in subNR_table['target'].unique():
			subname=subNR_table.loc[subNR_table['target'].eq(i)].iloc[0]
			Name=subname['taxname']+":"+subname['target']+"|"+str(subname['order'])
			try:
				handle=Entrez.efetch(db="protein", id=i, rettype="fasta_cds_na", retmode="text")
				record = SeqIO.read(handle, "fasta")
				print('>',re.sub(" ","_",Name),sep="",file=output)
				print(str(record.seq.translate()),file=output)
			except:
				print(i, " failed")





#Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";",index=False)
Output_table2=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results.tab",sep=";")


# Manually remove loci not from Filamentous thanks to the phylogeny checked by hand
nonfilamentous_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Table_with_nonfilamentous_loci.txt",sep="\t")

Output_table2['full_name']=Output_table2['species']+":"+Output_table2['target']+":"+Output_table2['tstart'].astype(int).astype(str)

Output_table2.loc[Output_table2['target'].str.contains("DWUI010000029.1"),"full_name"]="nan:DWUI010000029.1:3546"
Output_table2.loc[Output_table2['target'].str.contains("DWUI010000269"),"full_name"]="nan:DWUI010000269.1:803"

Output_table2=Output_table2.loc[~Output_table2['full_name'].isin(nonfilamentous_table['full_name'])]

#####################################################
# Find ORFs along the scaffolds with candidate loci #
#####################################################

# Write the scaffolds containing the hits into a file for each assembly 
tot=len(Output_table2['Assembly'].unique())
n=0
with open ("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.dna","a") as output:
	for assembly in Output_table2['Assembly'].unique():
		subOutput_table2=Output_table2.loc[Output_table2['Assembly'].eq(assembly)]
		assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes/"+assembly+".fna", "fasta"))
		for scaffold in subOutput_table2['target'].unique():
			#print(">",assembly,":",assembly_genome[scaffold].id,sep="",file=output)
			assembly_genome[scaffold].id=assembly+":"+assembly_genome[scaffold].id
			assembly_genome[scaffold].description=assembly+":"+assembly_genome[scaffold].description
			#print(textwrap.fill(str(assembly_genome[scaffold].seq),width=80),file=output)
			SeqIO.write(assembly_genome[scaffold],output,"fasta")
		n+=1
		print(n,"/",tot)

# Write the hits 
subprocess.run("cat /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Clusters/*.aa >> /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa" , shell=True)

# Run python ORF finder 
import pathlib  

Candidate_scaffolds="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.dna"
Candidate_scaffolds_ORFs_bed="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.bed"
Candidate_scaffolds_ORFs_dna="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.fna"
Candidate_scaffolds_ORFs_aa="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_scaffolds.aa"
outdir="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/"

subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/orfipy " +Candidate_scaffolds+"  --ignore-case  --procs 5 --outdir "+outdir+" --bed "+Candidate_scaffolds_ORFs_bed+" --min 150 --start ATG --dna "+Candidate_scaffolds_ORFs_dna, shell=True)

# Extract and translate the ORFs
ORF_bed=pd.read_csv(Candidate_scaffolds_ORFs_bed,sep="\t",header=None)
ORF_bed.columns=['Scaffold_name','ORF_start','ORF_end','ORF_name','zero','ORF_strand']

ORF_bed['ORF_name2']=ORF_bed['ORF_name'].str.replace(";.*","")
ORF_bed['ORF_name2']=ORF_bed['ORF_name2'].str.replace("ID=","")

with open(Candidate_scaffolds_ORFs_aa,"w") as output:
        record_dict=SeqIO.to_dict(SeqIO.parse(Candidate_scaffolds_ORFs_dna,"fasta"))
        for species in ORF_bed['Scaffold_name'].unique():
                for index, row in ORF_bed.loc[ORF_bed['Scaffold_name'].str.contains(species)].iterrows():
                        print(">",row['ORF_name2'],';',row['ORF_start'],"-",row['ORF_end'],"(",row['ORF_strand'],")",sep="",file=output)
                        print(record_dict[row['ORF_name2']].seq.translate(),file=output)

# Run mmseqs between ORFs and previously selected candidates

Predicted_ORFs_db=outdir+"scaffold_orfipy_db"
Filtred_loci="/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa"
subprocess.run("sed -i 's@ @_@g'  /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa", shell=True)

Filtred_loci_db=outdir+"All_candidate_hits_db"

Filtred_loci_vs_Predicted_ORFs_result=outdir+"Filtred_loci_vs_Predicted_ORFs_result"
Filtred_loci_vs_Predicted_ORFs_temp=outdir+"Filtred_loci_vs_Predicted_ORFs_tpm"
Filtred_loci_vs_Predicted_ORFs_table=outdir+"Filtred_loci_vs_Predicted_ORFs_result.m8"

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Candidate_scaffolds_ORFs_aa + " "+ Predicted_ORFs_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Filtred_loci + " "+ Filtred_loci_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  Filtred_loci_db  +" "+ Predicted_ORFs_db + " "+  Filtred_loci_vs_Predicted_ORFs_result + " "+ Filtred_loci_vs_Predicted_ORFs_temp + " -e 0.000001 --threads 30", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_loci_db +" "+ Predicted_ORFs_db + " "+  Filtred_loci_vs_Predicted_ORFs_result + " "+ Filtred_loci_vs_Predicted_ORFs_table + " --threads 10", shell=True)

ORF_vs_EVE_table=pd.read_csv(Filtred_loci_vs_Predicted_ORFs_table,sep="\t",header=None)
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
                print("outside")
                print(row['tstart']," : ", row['tend'])
                print(row['ORF_start']," : ", row['ORF_end'])
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

is_overlapped = lambda x: x['ORF_start'] >= x['ORF_end'].shift(fill_value=-1)
Output_table2_ORFs['ORF_overlapp_group'] = Output_table2_ORFs.sort_values(['full_name', 'ORF_start', 'ORF_end']) \
                .groupby(['full_name'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Output_table2_ORFs.loc[Output_table2_ORFs['ORF_start'].isna(),"ORF_overlapp_group"] = np.nan

# Calculate ORF_perc_best_hit 
Output_table2_ORFs['Best_hit_ORF_perc']= ((Output_table2_ORFs['ORF_end']-Output_table2_ORFs['ORF_start']) /3 )/Output_table2_ORFs['qlen']

Output_table2_ORFs.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results_ORFs.tab",sep=";",index=False)
Output_table2_ORFs=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ALL_results_ORFs.tab",sep=";")



Output_table2_ORFs['full_name'] =Output_table2_ORFs['species'].str.replace(" ","_")+"-"+Output_table2_ORFs['target']+"_"+Output_table2_ORFs['tstart'].astype(str)+"-"+Output_table2_ORFs['tend'].astype('str')+"("+Output_table2_ORFs['strand']+")|"+Output_table2_ORFs['order']
nonfilamentous_table=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Table_with_nonfilamentous_loci.txt",sep="\t")
Output_table2_ORFs=Output_table2_ORFs.loc[~Output_table2_ORFs['full_name'].isin(nonfilamentous_table['full_name'])]



#Output_table2_ORFs=Output_table2_ORFs.loc[Output_table2_ORFs['diff_evalue'].gt(min_control_diff_evalue)]

Output_table2_ORFs_sum=Output_table2_ORFs.drop_duplicates(subset = ["Assembly","target"])
Output_table2_ORFs_sum['Cumulative_scaffold_length'] = Output_table2_ORFs_sum.groupby(['Assembly'])['tlen'].transform('sum')
del Output_table2_ORFs['Cumulative_scaffold_length']
Output_table2_ORFs=Output_table2_ORFs.merge(Output_table2_ORFs_sum[['Assembly','Cumulative_scaffold_length']],on="Assembly")

# Generate binary table for the heatmap plot 

#manualyy fix AC81 for clavipes 
Output_table2_ORFs.loc[Output_table2_ORFs['Gene_name'].str.contains("ac81") & Output_table2_ORFs['species'].str.contains("clavipes"),"Best_hit_ORF_perc"]=0.90
Binary_Output_table2=Output_table2_ORFs.copy()
Binary_Output_table2 = Binary_Output_table2.drop_duplicates(subset = ["Assembly", "target","Gene_name"])
#Binary_Output_table2 = Binary_Output_table2.join(pd.crosstab(Output_table2['Assembly'], Output_table2['Gene_name']).add_prefix(''), on='Assembly')
Binary_Output_table2 = Binary_Output_table2.join(pd.crosstab(Output_table2_ORFs['Assembly'], Output_table2_ORFs['Gene_name']).add_prefix(''), on='Assembly')
Binary_Output_table2 = Binary_Output_table2.drop_duplicates()
Binary_Output_table2=Binary_Output_table2[['Assembly','full_name','tstart','tend','genus','species','family','order','Best_hit_ORF_perc','Cumulative_scaffold_length','Nb_EVEs','38K', 'ATPase', 'Integrase2', 'LbFVorf102', 'LbFVorf20', 'LbFVorf23', 'LbFVorf5', 'LbFVorf87', 'LbFVorf92','LbFVorf94', 'LbFVorf99', 'PDDEXK', 'ac38', 'ac81', 'dnapol', 'helicase', 'lcat', 'lef4', 'lef8', 'lef9', 'p33', 'p74', 'pif1', 'pif2', 'pif3', 'pif5']]
Binary_Output_table2=Binary_Output_table2.drop_duplicates(subset = "Assembly")
Binary_Output_table2.sort_values(by=["order","family","genus"], inplace = True)

# Add 2 if the ORF is not long enough compared to the best hit to be potentialy functional 

cols=['38K', 'ATPase', 'Integrase2', 'LbFVorf102', 'LbFVorf20', 'LbFVorf23', 'LbFVorf5', 'LbFVorf87', 'LbFVorf92','LbFVorf94', 'LbFVorf99', 'PDDEXK', 'ac38', 'ac81', 'dnapol', 'helicase', 'lcat', 'lef4', 'lef8', 'lef9', 'p33', 'p74', 'pif1', 'pif2', 'pif3', 'pif5']

for assembly in Binary_Output_table2['Assembly'].unique():
	subOutput_table2_ORFs=Output_table2_ORFs.loc[Output_table2_ORFs['Assembly'].eq(assembly)]
	subBinary_Output_table2=Binary_Output_table2.loc[Binary_Output_table2['Assembly'].eq(assembly)]
	for gene in cols:
		try:
			ORF_perc=subOutput_table2_ORFs.loc[subOutput_table2_ORFs['Gene_name'].eq(gene)]['Best_hit_ORF_perc'].iloc[0]
			if  ORF_perc >= 0.70:
				Binary_Output_table2.loc[Binary_Output_table2['Assembly'].eq(assembly),gene]=1
			else:
				Binary_Output_table2.loc[Binary_Output_table2['Assembly'].eq(assembly),gene]=2
		except:
			continue


#print(subBinary_Output_table2[gene].iloc[0]
#mask = Binary_Output_table2[cols].ne(0) & Binary_Output_table2['Best_hit_ORF_perc'].gt(0.85).to_numpy()[:, None]
#Binary_Output_table2[cols] = np.where(mask, 1, 2) * Binary_Output_table2[cols].ne(0)

#save 
Binary_Output_table2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Binary_filamentous_heatmap_table_EVEs_FL3.txt",sep=";",index=False)



#Output_table2_ORFs['Pseudogenized']="no"
#Filtred_loci_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/ORF_analysis/All_candidate_hits.aa","fasta"))










# Table output 
#Output_table2=Output_table2.loc[ ~(Output_table2['Cumulative_scaffold_length'].lt(200000) & Output_table2['Nb_EVEs'].gt(13))]

print(Output_table2[['Assembly','species','family','order','Nb_EVEs','Cumulative_scaffold_length']].drop_duplicates(subset = "Assembly"))
print(len(Output_table2[['Assembly','species','family','order','Nb_EVEs']].drop_duplicates(subset = "Assembly")))

