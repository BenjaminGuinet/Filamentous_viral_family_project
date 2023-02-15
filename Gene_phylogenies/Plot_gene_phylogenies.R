


library(ggtree)
library(ggplot2)



# Write all the core gene phylogenies 


files <- list.files(path="/Users/bguinet/Desktop/Filamentous_paper/Core_filamentous/Cluster_phylogenies", pattern="_AA.ali.treefile", full.names=TRUE, recursive=FALSE)


Figure_number=0

for (cluster in files){
  
  tree_file<-cluster
  tree_name<-gsub(".*/","",tree_file)
  
  tree_name<-gsub("\\..*","",tree_name)
  tree<-read.tree(tree_file)
  

  nb_taxa<-length(tree$tip.label)
  
  if (nb_taxa < 70){
    taxa_size=7
  }else{
    taxa_size=4
  }
  
  if (nb_taxa < 30){
    pdf_size_heigth=15
  }else{
    pdf_size_heigth=20
  }
  lim_right=10
  
  p1<- ggtree(tree)
  
  if (nrow(p1$data[grepl('WSSV',p1$data$label),])>0){
    tree<-root(tree, p1$data$label[grepl('WSSV',p1$data$label)])
  }else{
    tree<-midpoint.root(tree)
  }
  p1<- ggtree(tree)
  
  # By default let assume the loci s from a new filamentous virus 
  p1$data$category<-NA
  
  p1$data$category[p1$data$isTip=="TRUE"] <- "Non-Filamentous_Virus"
  p1$data$category[grepl('EfFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('CcFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('PoFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('PcFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('LbFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('LhFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('DFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('Dolicho',p1$data$label)]<-"Filamentous_Virus"
  
  
  Phylogeny_name<-tree_name
  Figure_number<- Figure_number+1
  
  mycolors=c("Filamentous_Virus"="#035397",
             "Non-Filamentous_New_Virus"="#FAAD80",
             "New_Filamentous_Virus"="#39A2DB",
             "Non-Filamentous_Virus"="#BD4B4B")
  
  
  p2 <- ggtree(p1$data, layout='rectangular') +
    geom_tippoint(
      mapping = aes(color = category),          # tip color by phyla. 
      size = 4,
      show.legend = FALSE) +
    geom_nodelab(size = 5,nudge_x = 0.1)  +
    geom_tiplab(    mapping = aes(color = category),                      # adds name of phyla to tip of its branch
                    offset = 0.5,
                    size = taxa_size,
                    #geom = "label",
                    align = TRUE,
                    face = "bold",
                    #label.size = 5,
                    #label.padding = unit(padding, "lines"), # amount of padding around the labels
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
  ggsave(paste0("/Users/bguinet/Desktop/Filamentous_paper/Core_filamentous/Phylogenies_pdf ",tree_name,".pdf"), width = 20, height = pdf_size_heigth, units = "in", limitsize = FALSE)
}


#############################################
# Write all the phylogenies with the insect EVE loci 

files <- list.files(path="/Users/bguinet/Desktop/Filamentous_paper/Core_filamentous/ALL_trees", pattern=".ali.treefile", full.names=TRUE, recursive=FALSE)

non_filamentous<-read.table("/Users/bguinet/Desktop/Filamentous_paper/Endogenous_FL_filamentous/Table_with_nonfilamentous_loci.txt",sep="\t",h=T)

Figure_number=0

for (cluster in files){
  
  tree_file<-cluster
  tree_name<-gsub(".*/","",tree_file)
  
  tree_name<-gsub("\\..*","",tree_name)
  tree<-read.tree(tree_file)
  
  #tree<-drop.tip(tree, tree$tip.label[grepl ("filamentous", tree$tip.label) ])
  
  nb_taxa<-length(tree$tip.label)
  
  if (nb_taxa < 70){
    taxa_size=7
  }else{
    taxa_size=4
  }
  
  if (nb_taxa < 30){
    pdf_size_heigth=15
  }else{
    pdf_size_heigth=20
  }
  lim_right=10
  
  p1<- ggtree(tree)
  
  if (nrow(p1$data[grepl('WSSV',p1$data$label),])>0){
    tree<-root(tree, p1$data$label[grepl('WSSV',p1$data$label)])
  }else{
  tree<-midpoint.root(tree)
  }
  p1<- ggtree(tree)
  
  # By default let assume the loci s from a new filamentous virus 
  p1$data$category<-NA
  
  p1$data$category[p1$data$isTip=="TRUE"] <- "Non-Filamentous_Virus"
  p1$data$category[grepl('\\|',p1$data$label)]<-"New_Filamentous_Virus"
  p1$data$category[grepl('_-_',p1$data$label)]<-"New_Filamentous_Virus"
  p1$data$category[grepl('___',p1$data$label)]<-"New_Filamentous_Virus"
  
  p1$data$label<-gsub("___","(+)",p1$data$label)
  p1$data$label<-gsub("_-_","(-)",p1$data$label)
  
  # Add non-filamentous for loci not found within a filamentous clade 
  subnon_filamentous<-non_filamentous[non_filamentous$Gene_name==tree_name,]
  
  for (loci in subnon_filamentous$full_name){
    p1$data$category[p1$data$label==loci]<-"Non-Filamentous_New_Virus"
  }
  
  #p1$data$category[is.na(p1$data$category)]<-"Virus"
  p1$data$category[grepl('EfFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('CcFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('PoFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('PcFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('LbFV_',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('LhFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('DFV',p1$data$label)]<-"Filamentous_Virus"
  p1$data$category[grepl('Dolicho',p1$data$label)]<-"Filamentous_Virus"
  
  
  Phylogeny_name<-tree_name
  Figure_number<- Figure_number+1
    
  mycolors=c("Filamentous_Virus"="#035397",
             "Non-Filamentous_New_Virus"="#FAAD80",
             "New_Filamentous_Virus"="#39A2DB",
             "Non-Filamentous_Virus"="#BD4B4B")
    
  
    p2 <- ggtree(p1$data, layout='rectangular') +
      geom_tippoint(
        mapping = aes(color = category),          # tip color by phyla. 
        size = 4,
        show.legend = FALSE) +
      geom_nodelab(size = 5,nudge_x = 0.1)  +
      geom_tiplab(    mapping = aes(color = category),                      # adds name of phyla to tip of its branch
        offset = 0.5,
        size = taxa_size,
        #geom = "label",
        align = TRUE,
        face = "bold",
        #label.size = 5,
        #label.padding = unit(padding, "lines"), # amount of padding around the labels
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
    ggsave(paste0("/Users/bguinet/Desktop/Filamentous_paper/Endogenous_FL_filamentous/Phylogenies_pdf ",tree_name,".pdf"), width = 20, height = pdf_size_heigth, units = "in", limitsize = FALSE)
}

