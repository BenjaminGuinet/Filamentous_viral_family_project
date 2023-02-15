import pandas as pd 
#tab['Lifestyle'] = np.nan
tab=pd.read_csv("All_taxonomy_informations.txt",sep=";")

tab_bippa=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Genomes_BIPAA_LBBE/Species_LBBE_BIPPA_list.txt")
tab_bippa['order']="Hymenoptera"
tab=pd.concat([tab,tab_bippa])
info_tab=pd.read_csv("/beegfs/data/bguinet/Filamentoviridae_vs_ALL_genomes/Bguinet_species_informations2.txt",sep="\t")



tab['species']=tab['species'].str.replace(" ","_")

tab['Lifestyle']=np.nan
for index, row in tab.iterrows():
  #print(row['species'])
  try:
    if len(info_tab.loc[info_tab['species'].str.contains(row['species'])]['Lifestyle'])>0:
      lifestyle=info_tab.loc[info_tab['species'].str.contains(row['species'])]['Lifestyle'].iloc[0]
      tab.loc[tab['species'].eq(row['species']),"Lifestyle"]=lifestyle
  except:
    continue 
  try:
      if len(info_tab.loc[info_tab['species'].str.contains(row['Accession'])]['Lifestyle'])>0:
        lifestyle=info_tab.loc[info_tab['species'].str.contains(row['Accession'])]['Lifestyle'].iloc[0]
        tab.loc[tab['Accession'].eq(row['Accession']),"Lifestyle"]=lifestyle
  except:
      continue 


tab.loc[tab['superfamily'].str.contains("Apoidea|Vespoidea|Formicoidea",na=False),"Lifestyle"]="freeliving"
#tab.loc[tab['superfamily'].str.contains("Apoidea|Vespoidea|Formicoidea",na=False)]
#tab.loc[~tab['superfamily'].str.contains("Apoidea|Vespoidea|Formicoidea",na=False) & tab['order'].str.contains("Hymeno")]

Number_free = len(tab.loc[tab['order'].str.contains("Hymeno") & tab['Lifestyle'].str.contains("free")]['species'].unique())
Number_endo = len(tab.loc[tab['order'].str.contains("Hymeno") & tab['Lifestyle'].str.contains("endo")]['species'].unique())
Number_ecto = len(tab.loc[tab['order'].str.contains("Hymeno") & tab['Lifestyle'].str.contains("ecto")]['species'].unique())

Number_free 
Number_endo
Number_ecto
