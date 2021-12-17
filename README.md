# Filamentous_viral_family_project

This project is a collaboration between the IMPI team in Tours (France) and the GEI team in Lyon (France). 
The objective of this project is to caracterize a new family of filamentous viruses.


Nous avons en notre possession 9 Génomes de virus filamenteux dont trois sont circularisés. L'idée de présent pipeline est de caractériser les gènes de ces génomes virales,
et de les annoter. 

## Recherche des ORFs dans les génomes 

** Programm : 
- getorf (version EMBOSS:6.6.0.0)
- orf_predictor.R 

Only open reading frames (ORFs) starting with a methionine and ending with a stop codon, with at least 50 amino acids and with minimal overlap (<23 nucleotides)
were considered as valuable candidates for being true ORFs
