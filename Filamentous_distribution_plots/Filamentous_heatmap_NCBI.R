


# Read the tree 
library(ape)
library(phytools)
library(phylogram)
library(ade4)
library(tidyverse)
library(ggplot2)


# Open the file with contains all the binary information of presence/absence of endogenous/exogenous FV among Diptera/Hymenoptera/Lepidoptera

Cluster <- read.table("/Users/bguinet/Desktop/Filamentous_paper/Binary_filamentous_heatmap_table_EVEs_FL.txt",sep=";",h=T)

#Cluster[!Cluster$species=="Encarsia formosa",]

Cluster <- Cluster %>%
  arrange(desc(species))%>%
  arrange(order) 

Cluster$pif1[Cluster$species=="Leptopilina boulardi"]<-0

Cluster_save<-Cluster 
Cluster$species[Cluster$species=="Junonia sp. 2 JZ-2019"]<-"Junonia sp"
Cluster$species[Cluster$species=="Dusona sp. PSUC_FEM 10030013"]<-"Dusona sp"
Cluster$species[Cluster$species=="Lissonota sp. PSUC_FEM 10030012"]<-"Lissonota sp"
Cluster$species[Cluster$species=="Dolichomitus sp. PSUC_FEM 10030005"]<-"Dolichomitus sp"
Cluster$species[Cluster$species=="Pseudomyrmex sp. PSW-54"]<-"Pseudomyrmex sp"


# put 1 if at least one assembly had the EVE
library(dplyr)

#cols<-c('X38K', 'ATPase', 'Integrase2', 'LbFVorf102', 'LbFVorf20', 'LbFVorf23', 'LbFVorf5', 'LbFVorf87', 'LbFVorf92', 'LbFVorf99', 'PDDEXK', 'ac38', 'ac81', 'dnapol', 'helicase', 'lcat', 'lef4','lef5', 'lef8', 'lef9', 'p33', 'p74', 'pif1', 'pif2', 'pif3', 'pif5','LbFVorf94','LbFVorf54')

cols<-c('LbFVorf23','LbFVorf54','LbFVorf87','LbFVorf94','LbFVorf99','LbFVorf5','LbFVorf20','LbFVorf102','lef4','lef5', 'lef8', 'lef9','dnapol', 'helicase','Integrase2','PDDEXK','LbFVorf92','p33','ac81','X38K','p69','ac38','lcat','ATPase','p74', 'pif1', 'pif2', 'pif3', 'pif5')

#Cluster <- Cluster %>%
#  group_by(species) %>%
#  mutate(across(all_of(cols), ~ if(1 %in% .x) replace(.x, is.na(.x)|
#                                                     .x %in% 0, 1) else .x)) %>%
#  ungroup



#Replace cumulative scaff length by the highest length if multiple same species 

#Cluster <- Cluster %>% 
#  group_by(species) %>% 
#  mutate(Cumulative_scaffold_length=max(Cumulative_scaffold_length,na.rm=TRUE))

Cluster<-as.data.frame(Cluster)


#Cluster<-Cluster[Cluster$Assembly %in%   Cluster[!duplicated(Cluster$species),]$Assembly,]

Cluster$pass_nb_EVEs <- "#EB5F34" #likely free-living virus 
Cluster$pass_nb_EVEs[Cluster$Nb_EVEs> 13 & Cluster$Cumulative_scaffold_length> 300000 ]<-"#35C44E9C" # Likely endogenized
Cluster$pass_nb_EVEs[Cluster$Nb_EVEs< 13 &  Cluster$Cumulative_scaffold_length > 300000]<-"#35C44E9C" # Likely endogenized
Cluster$pass_nb_EVEs[Cluster$Nb_EVEs< 13  & Cluster$Cumulative_scaffold_length < 300000]<-"#C3E477"


Cluster$genus <- NULL
Cluster$family <- NULL
Cluster$order <- NULL

Cluster<-Cluster[!duplicated(Cluster$species),]

Cluster_size_pass<-Cluster$pass_nb_EVEs

Cluster$pass_nb_EVEs<-NULL
Cluster$Cumulative_scaffold_length <- NULL
Cluster$Nb_EVEs <- NULL
Cluster$Best_hit_ORF_perc <- NULL
Cluster$Assembly <- NULL

rownames(Cluster) <- Cluster$species
Cluster$species <- NULL
Cluster$full_name<-NULL
Cluster$tstart<-NULL
Cluster$tend<-NULL

#Replace 
#only complete = 1  (1) 
#only incomplete = 2  (4)
#only pseudo = 3 (2 | 3 | 2.3 | 2.4 | 3.4)
#incomplete + complete = 4 (1.4)
#pseudo + complete = 5 (1.2 | 1.3 | 1.2.3 | 1.2.4 | 1.3.4)

# pseudo (2 et 3)
# incomplete without stop (4)
# complete (1)


# 6 means only pseudo 
Cluster[Cluster =="3"]<-6
Cluster[Cluster =="2"]<-6
Cluster[Cluster =="2-3"]<-6

# 5 means pseudo and complete
Cluster[Cluster =="1-2"]<-5
Cluster[Cluster =="1-3"]<-5
Cluster[Cluster =="1-2-3"]<-5
Cluster[Cluster =="1-2-4"]<-5
Cluster[Cluster =="1-3-4"]<-5


# 3 means only  incomplete  with pseudo 
Cluster[Cluster =="2-4"]<-3
Cluster[Cluster =="3-4"]<-3

# 2 means only incomplete without pseudo
Cluster[Cluster =="4"]<-2

# 4 means incomplete and complete 
Cluster[Cluster =="1-4"]<-4

# 1 means only complete 
Cluster[Cluster =="1"]<-1


Cluster$count_0<-NULL
Cluster$count_4<-NULL
Cluster$count_1.2<-NULL
Cluster$count_1.2.3<-NULL
Cluster$count_1.3<-NULL
Cluster$count_1.3.4<-NULL
Cluster$count_1<-NULL
Cluster$count_1.2.4<-NULL
Cluster$count_1.4<-NULL
Cluster$count_2<-NULL
Cluster$count_2.3<-NULL
Cluster$count_3<-NULL
Cluster$count_3.4<-NULL

library(gplots)
par(mar=c(7,4,4,2)+0.1) 


#15-10
dev.off()

rownames=rownames(Cluster)
Cluster<-sapply(Cluster, as.numeric)
rownames(Cluster)<-rownames
# Add number of copies within cell 
data.frame(colSums(Cluster != 0))
pdf(file = "/Users/bguinet/Desktop/Filamentous_paper/Heatmap_NCBI4.pdf",   
    width = 20, # The width of the plot in inches
    height = 20) 


black > 1
grey > 2
#red > 3
violet > 4 
#810600 > 5 


# 1 means only complete [black]
# 2 means only incomplete without pseudo [grey]
# 3 means only incomplete or complete with pseudo  [grey|red]
# 4 means incomplete and complete [noir|grey]
# 5 means pseudo and complete [noir|red]
# 6 means only pseudo 

heatmap.2(as.matrix(Cluster[ , cols]),# main = 'Heatmap of all Positive control candidats to a viral domestication',
          #reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
          # order by branch mean so the deepest color is at the top
          #Colv =hc, #color the branch labels 
          #cellnote=Cluster_cell_note,
          #notecol="black",
          #notecex=1,
          dendrogram='none',   
          labRow=as.expression(lapply(rownames(as.matrix(Cluster)), function(a) bquote(italic(.(a))))),
          RowSideColors=Cluster_size_pass,
          Rowv=FALSE,
          Colv=FALSE,  #   "yellow"                  "green" "yellow
          col = c("white","#282828","#9a9a9a","#D4B4B4","purple","#810600","#b51515"),         # color pattern of the heatmap
          #colCol = Cluster_names_vector_color,   #The cells color depending on the variable chosen
          trace="none",              # hide trace
          density.info="none",       # hide histogram
          margins = c(10,15),         # margin on top(bottom) and left(right) side.
          cexRow=1.5, cexCol = 1.5,      # size of row / column labels
          lhei=c(1, 12),           ## plot layout heights (H) of the rows in the plot.
          lwid=c(1,2),            ## plot layout width (L) of the rows in the plot.
          colsep=c(1:10000),       #Where to add column seperation
          rowsep = c(1:100000),      #here to add row seperation
          sepwidth=c(0.001,0.001),
          sepcolor="#D8D8D8",      #The color to separate de cells
          xlab = "Cluster names",
          key = F,
          offsetCol = -0.1,
          offsetRow = -0.1,
          srtCol = 45,
          # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
)
dev.off()

# T


# BArplot




very_likely_diptera<-length(unique(Cluster_save$species[Cluster_save$Nb_EVEs< 13  & Cluster_save$Cumulative_scaffold_length > 300000 & Cluster_save$order=="Diptera"]))
very_likely_hymeno<-length(unique(Cluster_save$species[Cluster_save$Nb_EVEs< 13  & Cluster_save$Cumulative_scaffold_length > 300000 & Cluster_save$order=="Hymenoptera"]))+1
very_likely_lepido<-length(unique(Cluster_save$species[Cluster_save$Nb_EVEs< 13  & Cluster_save$Cumulative_scaffold_length > 300000 & Cluster_save$order=="Lepidoptera"]))

likely_diptera<-length(unique(Cluster_save$species[Cluster_save$Nb_EVEs< 13  & Cluster_save$Cumulative_scaffold_length < 300000 & Cluster_save$order=="Diptera"]))
likely_hymeno<-length(unique(Cluster_save$species[Cluster_save$Nb_EVEs< 13  & Cluster_save$Cumulative_scaffold_length < 300000 & Cluster_save$order=="Hymenoptera"]))
likely_lepido<-length(unique(Cluster_save$species[Cluster_save$Nb_EVEs< 13  & Cluster_save$Cumulative_scaffold_length < 300000 & Cluster_save$order=="Lepidoptera"]))


# Correct by the number of genome input 

Nb_species=length(unique(Cluster_save$species))

#NCBI database species distrubtion
#Lepidoptera    892
#Diptera        358
#Hymenoptera    325

Diptera_expected<-   Nb_species*(358/(892+358+325))
Hymenoptera_expected<- Nb_species*(325/(892+358+325))
Lepidoptera_expected <- Nb_species*(892/(892+358+325))


dat <- data.frame(
  Status = c("very_likely","very_likely","very_likely","likely","likely","likely","Endogenized","Endogenized","Endogenized"), 
  Order=c("Diptera","Hymenoptera","Lepidoptera","Diptera","Hymenoptera","Lepidoptera","Hymenoptera","Diptera","Lepidoptera"), 
  Expected=c(Diptera_expected,Hymenoptera_expected,Lepidoptera_expected,Diptera_expected,Hymenoptera_expected,Lepidoptera_expected,0,0,0),
  Nb_species = c(very_likely_diptera,very_likely_hymeno,very_likely_lepido,likely_diptera,likely_hymeno,likely_lepido,6,0,0))


dat$perc_species[dat$Order=="Diptera"]<-    dat$Nb_species[dat$Order=="Diptera"]/(358)
dat$perc_species[dat$Order=="Hymenoptera"]<-    dat$Nb_species[dat$Order=="Hymenoptera"]/(325)
dat$perc_species[dat$Order=="Lepidoptera"]<-    dat$Nb_species[dat$Order=="Lepidoptera"]/(892)


ggplot(data=dat, aes(x=Order, y=perc_species, fill=Status)) +
  geom_bar(stat="identity",colour="black")+
  theme_bw()+ theme(axis.text.y = element_text(color="black", 
                                               size=7))+
  labs(x ="Species names", y = "Nb domesticated EVEs (dEVEs)")+
  #geom_point(size = 4.5,aes(y=Expected, shape="Exposure",fill="black"),shape=3,stroke = 1)+
  scale_fill_manual(values=c("#00335C","#C3E477", "35C44E9C","#056BF1"))+
  scale_y_continuous(n.breaks = 10)


ggsave(paste0("/Users/bguinet/Desktop/Filamentous_paper/Species_EVE_distribution_NCBI.pdf"), width = 5, height = 7, units = "in", limitsize = FALSE)


Count_diptera<-very_likely_diptera+likely_diptera
Count_hymeno<-very_likely_hymeno+likely_hymeno
Count_lepido<-very_likely_lepido+likely_lepido


M2 <- as.table(rbind(c(Hymenoptera_expected,Diptera_expected,Lepidoptera_expected), c(Count_hymeno, Count_diptera, Count_lepido))) # création d'une table 2 lignes/3 colonnes
dimnames(M2) <- list(gender=c("Nb species","Nb EVEs"), Genome_structure=c("Hymeno","Dipter", "Lepido"))# 

print(M2)
(test <- chisq.test(M2))
chisq.test(M2)$resi


Ecto_expected<-29
Endo_expected<-87
Free_expected<-209

Count_ecto<-0
Count_endo<-36
Count_free<-4
M3 <- as.table(rbind(c(Ecto_expected,Endo_expected,Free_expected), c(Count_ecto, Count_endo, Count_free))) # création d'une table 2 lignes/3 colonnes
print(M3)
(test <- chisq.test(M3))
chisq.test(M3)$resi

# Stats on most frequently observed core gene among EVEs 


Cluster_save<-Cluster_save[!duplicated(Cluster_save), ]
Cluster_save2<-t(Cluster_save)
Cluster_save2[Cluster_save2 != 0] = 1
Cluster_save2<-as.data.frame(Cluster_save2)
Cluster_save2[] <- sapply(Cluster_save2, as.numeric)

rowSums( Cluster_save2[,1:ncol(Cluster_save2)] )-1

Cluster_save2<-as.data.frame(Cluster_save2)%>% 
  mutate(sum = rowSums(across(where(is.numeric))))






