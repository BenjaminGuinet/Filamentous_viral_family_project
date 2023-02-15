


#This script will create a heatmap of all gene content between dsDNA genes 


# Read the tree 
library(ape)
library(phytools)
library(phylogram)

dsDNA_tree<-read.tree("/Users/bguinet/Desktop/Filamentous_paper/dsDNA_all-genes_phylogeny")


dsDNA_tree<-drop.tip(as.phylo(dsDNA_tree),"DFV")
dsDNA_tree<-drop.tip(as.phylo(dsDNA_tree),"DmSGHV")

dsDNA_tree$tip.label[dsDNA_tree$tip.label=="CpV"]<-"CpGV"
dsDNA_tree$tip.label[dsDNA_tree$tip.label=="tom"]<-"DmNV_tom"
dsDNA_tree$tip.label[dsDNA_tree$tip.label=="kal"]<-"DmNV_kal"
dsDNA_tree$tip.label[dsDNA_tree$tip.label=="mau"]<-"DmNV_mau"
dsDNA_tree$tip.label[dsDNA_tree$tip.label=="esp"]<-"DmNV_esp"

# Read the clusters table  

Cluster <- read.table("/Users/bguinet/Desktop/Filamentous_paper/All_core_gene_table.csv",sep=";",h=T,skip=1)
Cluster<-Cluster[!Cluster$X=="DmSGHV",]
Cluster <-subset(Cluster, select = -c(Odv.e66) )

names(Cluster) <- gsub(x = names(Cluster), pattern = "\\.", replacement = "-")  


Filamentous_unique<-c("LbFVorf23", "LbFVorf87", "LbFVorf94", "LbFVorf99")
Lefa_core<-c("Lef-4", "Lef-5", "Lef-8", "Lef-9" ,"helicase","Ac81")
Nadla_core<-c("pif-0-P74", "pif-1", "pif-2", "pif-3", "pif-5-odv-e56","DNApol","p33")
Hytro_fila<-c("LbFVorf5","LbFVorf20","LbFVorf54","LbFVorf102","lcat","LbFVorf92","PD--D-E-XK-nuclease")

Cluster<-Cluster %>% mutate_at(Filamentous_unique, funs(ifelse(.> 0, 3, .)))
Cluster<-Cluster %>% mutate_at(Lefa_core, funs(ifelse(.> 0, 4, .)))
Cluster<-Cluster %>% mutate_at(Nadla_core, funs(ifelse(.> 0, 5, .)))
Cluster<-Cluster %>% mutate_at(Hytro_fila, funs(ifelse(.> 0, 6, .)))



# Manually add sequences 
library(tidyverse)

List_names<-Cluster$X
Cluster$X<-NULL
Cluster<-as.matrix(Cluster)

rownames(Cluster)<-List_names

a=read.dendrogram(file = "/Users/bguinet/Desktop/Filamentous_paper/dsDNA_all-genes_phylogeny")

library(ape)
#Remove DFV 
a<-drop.tip(as.phylo(a),"DFV")
a<-drop.tip(as.phylo(a),"DmSGHV")
a<-as.dendrogram(a)
clade_order <- order.dendrogram(a)

clade_name <- labels(as.dendrogram(a))

clade_name[clade_name=="CpV"]<-"CpGV"
clade_name[clade_name=="tom"]<-"DmNV_tom"
clade_name[clade_name=="kal"]<-"DmNV_kal"
clade_name[clade_name=="mau"]<-"DmNV_mau"
clade_name[clade_name=="esp"]<-"DmNV_esp"

clade_position <- data.frame(clade_name, clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name,row.names(Cluster))
Cluster <- Cluster[new_order,]
clade_name <- labels(as.dendrogram(a))

clade_name[clade_name=="CpV"]<-"CpGV"
clade_name[clade_name=="tom"]<-"DmNV_tom"
clade_name[clade_name=="kal"]<-"DmNV_kal"
clade_name[clade_name=="mau"]<-"DmNV_mau"
clade_name[clade_name=="esp"]<-"DmNV_esp"
rownames(Cluster) <- clade_name

dev.off()

heatmap.2(as.matrix(Cluster), main = 'Heatmap of all Positive control candidats to a viral domestication',
          #reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
          # order by branch mean so the deepest color is at the top
          Rowv = as.dendrogram(as.phylo(dsDNA_tree)),
          #Colv =hc, #color the branch labels 
          notecol="white",
          notecex=2,
          Colv=FALSE,  #   "yellow"                  "green" "yellow
          col = c("white",
                  "#393E46",     
                  "darkgreen",
                  "#00ABB3",
                  "#6D9886",
                  "#FEDB39", "#FCF6C4"),         # color pattern of the heatmap
          #colCol = Cluster_names_vector_color,   #The cells color depending on the variable chosen
          trace="none",              # hide trace
          density.info="none",       # hide histogram
          margins = c(10,18),         # margin on top(bottom) and left(right) side.
          cexRow=0.8, cexCol = 0.8,      # size of row / column labels
          lhei=c(1, 12),           ## plot layout heights (H) of the rows in the plot.
          lwid=c(1,4),            ## plot layout width (L) of the rows in the plot.
          colsep=c(1:10000),       #Where to add column seperation
          rowsep = c(1:100000),      #here to add row seperation
          sepwidth=c(0.001,0.001),
          sepcolor="#C8C6C6",      #The color to separate de cells
          xlab = "Cluster names",
          key = F,
          offsetCol = -0.1,
          offsetRow = -0.1,
          srtCol = 45,
          # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
)



#15-10
dev.off()



### On all genes 
library(tidyr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ape)
library(phytools)
library(phylogram)
library(gplots)


dat <- read.table("/Users/bguinet/Desktop/Filamentous_paper/fichier_reduit_homologue_benjamin.csv",sep=";",h=T)

#dat <- dat[! is.na(dat$CcFV1),]
#dat <- dat[! is.na(dat$CcFV2),]
#dat<-dat[c('Groups','CcFV2','CcFV1')]

dat[] <- lapply(dat, as.character)

dat[dat == ''] <- NA

dat2<-dat %>%
  pivot_longer(-Groups) %>%
  group_by(Groups, Groups2 = factor(gsub("_.*", "", name))) %>%
  filter(!is.na(value)) %>%
  summarise(value = paste(paste(name, value, sep = "-"), collapse = ";")) %>%
  ungroup() %>%
  pivot_wider(names_from = "Groups",
              values_from = "value") %>%
  complete(Groups2)

#Merge column if duplicates in rows between columns
row_names<-dat2$Groups2
dat2$Groups2<-NULL

dat2<- dat2%>%
  rowid_to_column() %>%
  pivot_longer(-rowid) %>%
  filter(!is.na(value)) %>%
  group_by(rowid, value) %>%
  mutate(new_name = paste(name, collapse = "|")) %>%
  separate_rows(new_name, sep = "\\|") %>%
  group_by(name) %>%
  mutate(new_name = paste(unique(new_name), collapse = "|")) %>%
  group_by(value) %>%
  filter(nchar(new_name) == max(nchar(new_name))) %>%
  ungroup() %>%
  select(-name) %>%
  pivot_wider(names_from = new_name, values_from = value, values_fn = ~ paste(unique(.x), collapse = "|")) %>%
  complete(rowid = full_seq(c(1, rowid), 1))


dat2$rowid<-NULL
col_names <- paste("G", 0:length(dat2), sep = "")
colnames(dat2)<-col_names

dat3<-dat2

write.table(dat2, file="/Users/bguinet/Desktop/Filamentous_paper/All_gene_content_group_names.csv",sep=";")


# Replace non Nan by 1 


dat2<-replace(data.frame(lapply(dat2, as.character), stringsAsFactors = FALSE),
              !is.na(dat2), 1)

dat2<-replace(data.frame(lapply(dat2, as.character), stringsAsFactors = FALSE),
              is.na(dat2),0)
rownames(dat2)<-row_names

a=read.dendrogram(file = "/Users/bguinet/Desktop/Filamentous_paper/Filamentous_plus_Hytro.nwk")

library(ape)
#Remove DFV 
a<-as.dendrogram(a)
clade_order <- order.dendrogram(a)
clade_name <- labels(as.dendrogram(a))


clade_position <- data.frame(clade_name, clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name,row.names(dat2))

dat2 <- dat2[clade_name, ]  

clade_name <- labels(as.dendrogram(a))

dat2<-sapply(dat2, as.numeric)
rownames(dat2) <- clade_name



heatmap.2(as.matrix(dat2),
          reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
          # order by branch mean so the deepest color is at the top
          Rowv = as.dendrogram(as.phylo(a)),
          
          #Rowv = c("CcFV1",  "CcFV2",  "EfFV","LbFV",   "LhFV",   "PcFV",   "PoFV",  "MdSGHV","GpSGHV"),
          #Colv =hc, #color the branch labels 
          notecol="white",
          notecex=2,
          #Colv=FALSE,  #   "yellow"                  "green" "yellow
          col = c("white","#548CA8"),         # color pattern of the heatmap
          #colCol = Cluster_names_vector_color,   #The cells color depending on the variable chosen
          trace="none",              # hide trace
          density.info="none",       # hide histogram
          margins = c(10,18),         # margin on top(bottom) and left(right) side.
          cexRow=4, cexCol = 0.2,      # size of row / column labels
          lhei=c(1, 12),           ## plot layout heights (H) of the rows in the plot.
          lwid=c(1,4),            ## plot layout width (L) of the rows in the plot.
          colsep=c(1:10000),       #Where to add column seperation
          rowsep = c(1:100000),      #here to add row seperation
          sepwidth=c(0.0001,0.0001),
          sepcolor="#D8D8D8",      #The color to separate de cells
          xlab = "Cluster names",
          key = F,
          offsetCol = -0.1,
          offsetRow = -0.1,
          srtCol = 45,
          # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
)

50-20
