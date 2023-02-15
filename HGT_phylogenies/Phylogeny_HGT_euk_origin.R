# Midpoint root a tree 
library("phytools")
library(ggplot2)
library(ggtree)

files <- list.files(path="/Users/bguinet/Desktop/Filamentous_paper/Non_viral_phylogeny_cluster3/All_tree", pattern=".treefile", full.names=TRUE, recursive=FALSE)

#HGT_table<-read.table("/Users/bguinet/Desktop/Filamentous_paper/Non_viral_phylogeny_cluster2/HGT_ORFS_table.txt",sep="\t",h=T)

Group_name_table<-read.table("/Users/bguinet/Desktop/Filamentous_paper/Non_viral_phylogeny_cluster3/Group_name_matthieu_cluster_non_viral_correspondance.txt",sep=";",h=T)
Figure_number=0

for (cluster in files){
  
tree_file<-cluster
tree_name<-gsub(".*/","",tree_file)
Cluster_name<-gsub(".aa.ali.treefile","",tree_name)
Group_name<-Group_name_table$Group_name[Group_name_table$Cluster_hmmer==Cluster_name]
Function_name<-Group_name_table$Hmmer_function[Group_name_table$Cluster_hmmer==Cluster_name]
Phylogeny_name<-paste0(Group_name," (",Function_name,")" )

if (identical(Group_name, character(0)) | Cluster_name=="Cluster_JmJC"){
  print(Group_name)
  keep_it<-"no"
}else{
  Figure_number<- Figure_number+1
  keep_it<-"yes"
}

# Remove duplicated taxa 

tree_name<-gsub("\\..*","",tree_name)
tree<-read.tree(tree_file)

tree<-drop.tip(tree, tree$tip.label[grepl ("filamentous", tree$tip.label) ])

nb_taxa<-length(tree$tip.label)

if (nb_taxa < 150){
  taxa_size=4
 }else{
  taxa_size=2
 }

node_size=2
boot_siz=4
nudge = 0.07
if (nb_taxa > 200){
  taxa_size=1
  node_size=1
  boot_size=2
  nudge = 0.05
}

if (nb_taxa < 30){
  pdf_size_heigth=15
}else{
  pdf_size_heigth=20
}

if (grepl("Cluster389|Cluster324|Cluster279",tree_name)){
  lim_right=30
}else{
  lim_right=10
  }

tree<-midpoint.root(tree)

p1<- ggtree(tree)

p1$data$category <- NA

p1$data$category[grepl('Eukaryota',p1$data$label)]<-"Eukaryota"
p1$data$category[grepl('Bacter',p1$data$label)]<-"Bacteria"
p1$data$category[grepl('Arch',p1$data$label)]<-"Archea"
p1$data$category[is.na(p1$data$category)]<-"Virus"
p1$data$category[grepl('EfFV',p1$data$label)]<-"LbFV-like virus"
p1$data$category[grepl('CcFV',p1$data$label)]<-"LbFV-like virus"
p1$data$category[grepl('PoFV',p1$data$label)]<-"LbFV-like virus"
p1$data$category[grepl('PcFV',p1$data$label)]<-"LbFV-like virus"
p1$data$category[grepl('LbFV',p1$data$label)]<-"LbFV-like virus"
p1$data$category[grepl('LhFV',p1$data$label)]<-"LbFV-like virus"
p1$data$category[grepl('DFV',p1$data$label)]<-"LbFV-like virus"


# Continue only if filamentous viruses 
if (nrow(p1$data[p1$data$category=="LbFV-like virus",])>0){

  # Continue only if clear event 
  
  
  #library(rapport)
  #if (length(Function_name)>0){
  #  Phylogeny_name<-paste0(Cluster_name," (",Function_name,")" )
  #}else{
  #  Phylogeny_name<-paste0(Cluster_name)
  #}

mycolors=c("Eukaryota"="#54B435","Bacteria"="#FFB72B","Archea"="grey","Virus"="#D72323","LbFV-like virus"="#00ABB3")

p2 <- ggtree(tree, layout='rectangular',size=0.5) %<+% p1$data +
  geom_tippoint(
    mapping = aes(color = category),          # tip color by phyla. 
    size = node_size,
    show.legend = FALSE) +
  geom_nodelab(size = boot_size, col= "red",nudge_x = nudge)  +
  geom_tiplab(                          # adds name of phyla to tip of its branch
    aes(color = category),
    #color = 'black',                      
    offset = 0.5,
    fontface=2,
    size = taxa_size,
    #geom = "label",
    align = TRUE,
    face = "bold",
    #label.size = 0,
    #label.padding = unit(0.15, "lines"), # amount of padding around the labels
    linetype = "dashed") +
  ggtitle("Phylogenetic tree of Tert")+  # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) +
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors) + xlim(0,lim_right)  + ggtitle(paste0("FigureS",Figure_number," : ", Phylogeny_name))

                                                                      
print(Figure_number)                                                                
p2

if (keep_it=="yes"){
ggsave(paste0("/Users/bguinet/Desktop/Filamentous_paper/Non_viral_phylogeny_cluster3/Phylogenies_pdf/",paste0("FigureS",Figure_number,'_',Group_name,"_",Cluster_name),".pdf"), width = 20, height = pdf_size_heigth, units = "in", limitsize = FALSE)
}

if (keep_it=="no"){
  ggsave(paste0("/Users/bguinet/Desktop/Filamentous_paper/Non_viral_phylogeny_cluster3/Phylogenies_pdf/",Cluster_name,".pdf"), width = 20, height = pdf_size_heigth, units = "in", limitsize = FALSE)
}
  }
}

