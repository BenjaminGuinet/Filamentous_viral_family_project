

import subprocess
from taxadb.taxid import TaxID
from Bio import SeqIO 
import pandas as pd 
import os,re 
already_ran=[]

if os.path.exists("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/All_taxonomy_informations.txt"):
	Tax_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/All_taxonomy_informations.txt",sep=";")
	already_ran=Tax_tab['Assembly'].unique()
else:
	Tax_tab=pd.DataFrame(columns=['Taxid','Assembly', 'family','subfamily','species','superfamily','order','superorder'])


Accession_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/ALL_accession_numbers.txt",header=None)
Accession_tab.columns=['Accession']
Accession_tab=Accession_tab.loc[~Accession_tab['Accession'].isin(already_ran)]

Accession_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/ALL_accession_numbers.txt",header=None)
Accession_tab.columns=['Accession']

Accession_tab2=Accession_tab.merge(Tax_tab,left_on="Accession",right_on="Assembly",how="outer")
Accession_tab2=Accession_tab2.loc[~Accession_tab2['Accession'].isna()]


n_accessions=len(Accession_tab2.loc[Accession_tab2['Taxid'].isna()]['Accession'].unique())
n=0



#Function to return smaller lists of a larger list in order to question the database with several accession number instead of just one.
def chunk(items,n):
    for i in range(0,len(items),n):
        yield items[i:i+n]

number_of_sub_acc_number_lists=1

Acc_number_list=Accession_tab2.loc[Accession_tab2['Taxid'].isna()]['Accession'].unique()
# Break the larger list into 10 length lists.
count_row = len(Acc_number_list)
filecount = 1

def str2frame(estr, sep = '\t', lineterm = '\n', set_header = True):
    dat = [x.split(sep) for x in estr.split(lineterm)][1:-1]
    df = pd.DataFrame(dat)
    if set_header:
        df = df.T.set_index(0, drop = True).T # flip, set ix, flip back
    return df

taxids = TaxID(dbtype='sqlite', dbname='/beegfs/data/bguinet/taxadb2/taxadb.sqlite')
#For each sub-list of access numbers of size 10, retrieve the equivalent taxonomic IDs and stor it into taxids variable.
for sublist in chunk(Acc_number_list,number_of_sub_acc_number_lists):
	sub_accessions='|'.join(sublist)
	try:
        	taxid=subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/esearch -db assembly -q '"+ sub_accessions +"' | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid", shell=True,capture_output=True, text=True)
	        taxid=taxid.stdout
		taxid_tab=str2frame(taxid)
		#print(taxid_tab)
		taxid_tab=taxid_tab.T.reset_index().T.reset_index(drop=True)
		taxid_tab.columns=['Accession','Taxid']
		print(taxid_tab)
	except:
		print(sub_accessions)
	#print(taxid_tab)
	for index, row in taxid_tab.iterrows():
		try:
			Taxid_dictionnary = {}
			Taxid_dictionnary[row['Taxid']]= taxids.lineage_name(row['Taxid'],ranks=True,reverse=True)
			taxid_db=pd.DataFrame({k: dict(v) for k, v in Taxid_dictionnary.items()}).T
			taxid_db['Taxid']=row['Taxid']
			taxid_db['Assembly']=row['Accession']
			print(taxid_db)
			Tax_tab=Tax_tab.append(taxid_db)
		except:
			taxid_db['Taxid']=row['Taxid']
			taxid_db['Assembly']=row['Accession']
			taxid_db['species']="did_not_work"
			Tax_tab=Tax_tab.append(taxid_db)
        print(filecount, " / ", count_row)
        filecount+=100
	Tax_tab.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/All_taxonomy_informations.txt",sep=";",index=False)


Accession_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_NCBI/ALL_accession_numbers.txt",header=None)
Accession_tab.columns=['Accession']

Accession_tab2=Accession_tab.merge(Tax_tab,left_on="Accession",right_on="Assembly",how="outer")

Accession_tab2=Accession_tab2.loc[~Accession_tab2['Accession'].isna()]


for accession in Accession_tab2.loc[Accession_tab2['Taxid'].isna()]['Accession']:
	try:
		taxid=subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/esearch -db assembly -q '"+ accession +"' | esummary | xtract -pattern DocumentSummary -element AssemblyAccession,Taxid", shell=True,capture_output=True, text=True)
		taxid=taxid.stdout
		taxid = re.sub(".*\t","",taxid)
		taxid = re.sub("\n.*","",taxid)
		Taxid_dictionnary = {}
		Taxid_dictionnary[taxid]= taxids.lineage_name(taxid,ranks=True,reverse=True)
		taxid_db=pd.DataFrame({k: dict(v) for k, v in Taxid_dictionnary.items()}).T
		taxid_db['Taxid']=taxid
		taxid_db['Assembly']=accession
		taxid_db
		Tax_tab=Tax_tab.append(taxid_db)
	except:
		print(accession)



Accession_tab2 = Accession_tab2.drop_duplicates()
Accession_tab2.to_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/All_taxonomy_informations.txt",sep=";",index=False)

#Count number order 
Accession_tab2 = Accession_tab2.drop_duplicates(subset = "species")
Accession_tab2.order.value_counts()
