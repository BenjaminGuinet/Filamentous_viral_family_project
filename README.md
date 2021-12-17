# Filamentous_viral_family_project

This project is a collaboration between the IMPI team in Tours (France) and the GEI team in Lyon (France). 
The objective of this project is to caracterize a new family of filamentous viruses.


We have in our possession 9 filamentous virus genomes, three of which are circularised. 
The idea of this pipeline is to :
- (i) characterize the genes of these viral genomes
- (ii) annotate them
- (iii) reconstruct the phylogeny of this group within the dsDNA virus clade

## Search for ORFs in the genomes 

**Programs** : 
- getorf (version EMBOSS:6.6.0.0)
- orf_predictor.R 

Only open reading frames (ORFs) starting with a methionine and ending with a stop codon, with at least 50 amino acids and with minimal overlap (<23 nucleotides)
were considered as valuable candidates for being true ORFs

**Output files**

- Virusname_prediction_option_1.pdf (Genomic plot of the filtred predicted ORFs for each contigs)
- Virusname_prediction_option_1.gff (Gff file with coordinates of each filtred ORFs)
- Virusname_prediction_option_1.fa (Fasta file with predicted filtred ORFs in AminoAcide format) 
- Virusname_getorf_option_1.fa (Fasta file with predicted ORFs from getorf)
 

## Clustering of homologous ORFs

