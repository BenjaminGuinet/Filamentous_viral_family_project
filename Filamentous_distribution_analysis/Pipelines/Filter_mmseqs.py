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
#python3 Filter_mmseqs.py -f /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/GCA_905404225.1_result.m8  -o /beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/GCA_905404225.1_result_filtred.m8

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
        Mmseqs_tab=pd.read_csv(file,sep="\t",header=None,encoding= 'ISO-8859-1')
        #if  filename in ['GCA_019393585.1_result_orf94.m8','GCA_015476485.1_result_orf94.m8','GCA_011634795.1_result_orf94.m8']:
        #   Mmseqs_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','tstart','tend','qstart','qend','evalue','bits','qlen','tlen','tcov','qcov'] 
        Mmseqs_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
        print(filename)
        print(Mmseqs_tab)
        #Filter 
        #Keep small qcov if multiple ORFs along a contig 
        #Mmseqs_tab=Mmseqs_tab.loc[Mmseqs_tab['qcov'].ge(0.25)]
        Mmseqs_tab=Mmseqs_tab.loc[Mmseqs_tab['bits'].ge(40)]
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
           print("did not work")
          Mmseqs_tab.to_csv(output,sep=";",index=False)
         else:
          Output_table=pd.DataFrame(columns=['Taxid', 'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'Assembly', 'new_start', 'new_end', 'strand', 'Taxid_x', 'clade', 'class', 'cohort', 'family', 'genus', 'infraclass', 'infraorder', 'kingdom', 'no rank', 'order', 'parvorder', 'phylum', 'species', 'subclass', 'subfamily', 'suborder', 'subphylum', 'superfamily', 'superkingdom', 'superorder', 'Taxid_y'])
          Output_table.to_csv(output,sep=";",index=False) 
      else:
          Output_table=pd.DataFrame(columns=['Taxid', 'query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'Assembly', 'new_start', 'new_end', 'strand', 'Taxid_x', 'clade', 'class', 'cohort', 'family', 'genus', 'infraclass', 'infraorder', 'kingdom', 'no rank', 'order', 'parvorder', 'phylum', 'species', 'subclass', 'subfamily', 'suborder', 'subphylum', 'superfamily', 'superkingdom', 'superorder', 'Taxid_y'])
          Output_table.to_csv(output,sep=";",index=False)
 
