#Compute patristic distance between sub groups within a phylogny

library(ape)
tree<-read.tree("/Users/bguinet/Desktop/Papier_scientifique/dsDNA_all_genes_iqtree.nwk")

distance_matrix<-cophenetic.phylo(tree)
patristic_tab<-as.data.frame(as.table(distance_matrix))
colnames(patristic_tab)<-c("C1","C2","Distance")

write.table(patristic_tab,file="/Users/bguinet/Desktop/Filamentous_paper/dsDNA_all_core_iqtree_patristic_distances.tab",sep=";")

patristic_tab<- patristic_tab[!patristic_tab$C1 == patristic_tab$C2,]
patristic_tab<-patristic_tab[!duplicated(patristic_tab$Distance),]


Baculoviridae<-c("LdMNPV","AcMNPV","CpV","NeseNPV","CuniNPV")

Nudiviridae<-c("GbNV","OrNV","tom","esp","kal","DiNV","mau","ToNV","HzNV1","HzNV2","DhNV","PmNV","HgNV")
Alpha_Nudiviridae<-c("GbNV","mau","DiNV","tom")
Beta_Nudiviridae<-c("HzNV1","HzNV2")
Delta_Nudiviridae<-c("ToNV")
Gamma_Nudiviridae<-c("DhNV","PmNV")


LbFV<-c("LbFV","LhFV","EfFV","PoFV","PcFV","CcFV2","CcFV1")

Hytrosaviridae<-c("GpSGHV","MdSGHV","DmSGHV")

All_sequences <-c("LdMNPV","AcMNPV","CpV","NeseNPV","CuniNPV","GbNV","OrNV","tom","esp","kal","DiNV","mau","ToNV","HzNV1","HzNV2","DhNV","PmNV","HgNV",
                   "LbFV","LhFV","EfFV","PoFV","PcFV","CcFV2","CcFV1","GpSGHV","MdSGHV","DmSGHV" )


library(tidyverse)
lists1 <- stack(lst(Baculoviridae,Nudiviridae,Hytrosaviridae,LbFV))
patristic_tab1<-bind_cols(patristic_tab, 
          patristic_tab %>% 
            map_dfc(~ lists1$ind[match(.x, lists1$values)]) %>% 
            unite(col = "Lists", sep = "|"))
lists2 <- stack(lst(Baculoviridae,Alpha_Nudiviridae,Beta_Nudiviridae,Delta_Nudiviridae,Gamma_Nudiviridae,Hytrosaviridae,LbFV))
patristic_tab2<-bind_cols(patristic_tab, 
                          patristic_tab %>% 
                            map_dfc(~ lists2$ind[match(.x, lists2$values)]) %>% 
                            unite(col = "Lists", sep = "|"))
#Remove if NA 
patristic_tab<-rbind(patristic_tab1,patristic_tab2)

patristic_tab<-patristic_tab[!grepl("NA\\|",patristic_tab$Lists),]
patristic_tab$Lists<-gsub("\\|NA","",patristic_tab$Lists)



patristic_tab<-patristic_tab[!patristic_tab$Distance==0,]

library(dplyr)
library(ggplot2)

# ALLL 
patristic_tab %>%
  mutate(Lists = fct_reorder(Lists, Distance, .fun='median')) %>%
  ggplot( aes(x=reorder(Lists, Distance), y=Distance, fill=Lists)) + 
  geom_boxplot() +
  theme_bw(base_size = 14) +
  guides(fill=guide_legend("family comparisons"))+
  xlab("") +
  ylab("Patristic distance") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
scale_fill_discrete(guide = guide_legend(title = "Families/Families comparisons"))+
  scale_fill_manual(values=c("Baculoviridae|Baculoviridae"="#ece30d", "Nudiviridae|Baculoviridae"="#ffa500", "Hytrosaviridae|Baculoviridae"="#008080", "Hytrosaviridae|Hytrosaviridae"="#8692bb", "Hytrosaviridae|Nudiviridae"="#16B336","Baculoviridae|LbFV"="#0061FE",
                             "Hytrosaviridae|LbFV"="#1671B3","LbFV|LbFV"="#69c6e0", "Nudiviridae|LbFV"="#808671","Nudiviridae|Nudiviridae"="#935215"))


# Only families

patristic_tab[patristic_tab$Lists %in% c("LbFV|LbFV","Baculoviridae|LbFV","Nudiviridae|LbFV" ,"Hytrosaviridae|LbFV" ,"Baculoviridae|Baculoviridae","Nudiviridae|Baculoviridae" ,"Hytrosaviridae|Baculoviridae","Nudiviridae|Nudiviridae","Hytrosaviridae|Nudiviridae","Hytrosaviridae|Hytrosaviridae"),] %>%
  mutate(Lists = fct_reorder(Lists, Distance, .fun='median')) %>%
  ggplot( aes(x=reorder(Lists, Distance), y=Distance, fill=Lists)) + 
  geom_boxplot() +
  theme_bw(base_size = 14) +
  guides(fill=guide_legend("family comparisons"))+
  xlab("") +
  ylab("Patristic distance") + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_discrete(guide = guide_legend(title = "Families/Families comparisons"))+
  scale_fill_manual(values=c("Baculoviridae|Baculoviridae"="#ece30d", "Nudiviridae|Baculoviridae"="#ffa500", "Hytrosaviridae|Baculoviridae"="#008080", "Hytrosaviridae|Hytrosaviridae"="#8692bb", "Hytrosaviridae|Nudiviridae"="#16B336","Baculoviridae|LbFV"="#0061FE",
                             "Hytrosaviridae|LbFV"="#1671B3","LbFV|LbFV"="#69c6e0", "Nudiviridae|LbFV"="#808671","Nudiviridae|Nudiviridae"="#935215"))

