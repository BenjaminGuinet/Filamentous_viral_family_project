import csv
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database() #if needed, to upload the taxonomy
data_dict = {}


with open('Filamentous_distribution_taxid_among_insects_nodup.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=";")
    next(reader) # skip first row
	for row in reader:
		key = row[1];
		data_dict[key]=row

taxids=list(data_dict.keys())
tree = ncbi.get_topology(taxids)
