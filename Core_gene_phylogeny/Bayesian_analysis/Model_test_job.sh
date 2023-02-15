
perl /beegfs/home/bguinet/these_scripts_2/catfasta2phyml.pl -f --concatenate --verbose *_AA.ali.trimmed > /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis2/Concatenated_sequences_dsDNA.aln  2>  /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis/partition_dsDNA.txt
grep '_AA.ali.trimmed =' /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis/partition_dsDNA.txt >> /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis2/partition_dsDNA.tab
sed -i "s@Cluster@AA, Cluster@g" /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis2/partition_dsDNA.tab

#Run Iqtree model finder

/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis2/Concatenated_sequences_dsDNA.aln -spp /beegfs/data/bguinet/LbFV_family_project/Clustering/Cluster_phylogeny_core/Mrbayes_analysis2/partition_dsDNA.tab -m TESTONLY -nt 10

