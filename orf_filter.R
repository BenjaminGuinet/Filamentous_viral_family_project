#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#  argument  needed :
# the genome.fa,
# minorf size (in nucleotides) used by getorf,
# maxoverlap (in nucl.) authorized ,
# detection method used by getorf (default 1 : should start with Met and end with a stop)

#Usage example :

# /beegfs/data/soft/R-4.0.5/bin/Rscript orf_predictor.R /beegfs/data/bguinet/LbFV_family_project/Genomes/LbFV_free/LbFV.fna 150 23 1
# load libraries
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(gggenes)
library(ggplot2)
library(gridExtra)

#######

# import data

# virus genome
genome=readDNAStringSet(args[1])
contig_length=width(genome)
names(contig_length)=strsplit(names(genome), " ")[[1]][1]
names(genome) <- strsplit(names(genome), " ")[[1]][1]

# getorf predictions using getorf from embo
output_getorf=paste0(tools::file_path_sans_ext(args[1]), paste0("_getorf_option", args[4],".fa"))
command=paste0("getorf -sequence ", args[1]," -minsize ", args[2], " -find ", args[4]," -outseq ", output_getorf)
print(command)
system(command = command)

# load getorf predictions:
orf_pred=readBStringSet(output_getorf)

# Extract orf positions from getorf fasta file
names=names(orf_pred)
#length(names)
names=strsplit(names, " ")



names=lapply(names, FUN=function(x) {
  contig=(x[1]) # extract contig name
  start=sub(pattern = "[", replacement = "", x = x[2], fixed = T)
  stop=sub(pattern = "]", replacement = "", x = x[4], fixed = T)
  res=c(contig, start, stop)
  return(res)
}
)

orf_coord=do.call(rbind.data.frame, names)
names(orf_coord)=c("contig", "start", "stop")
orf_coord$start=as.numeric(orf_coord$start)
orf_coord$stop=as.numeric(orf_coord$stop)
orf_coord$strand="+"
orf_coord$strand[apply(X = orf_coord, MARGIN = 1, FUN = function(x){x[2]>x[3]})]="-"
orf_coord$start2=as.numeric(apply(X = orf_coord, 1, FUN=function(x){min(x[2:3])}))
orf_coord$stop2=as.numeric(apply(X = orf_coord, 1, FUN=function(x){max(x[2:3])}))
dim(orf_coord)
orf_coord$contig<-gsub("_[^_]+$","",orf_coord$contig)
################################################@
# FILTER ORFs
################################################@


gr_orfs=GRanges(ranges = IRanges(orf_coord$start2,orf_coord$stop2), strand = orf_coord$strand, seqnames = orf_coord$contig, seqinfo = contig_length)
dim(gr_orfs)
# find the predicted orfs entirely nested within others one and remove them
within_orfs = findOverlaps(query = gr_orfs, type = "within", drop.self=TRUE,  ignore.strand =TRUE)
if (length(within_orfs)>0){ # if some are nested, remove them
gr_orfs_subset=gr_orfs[-queryHits(within_orfs)]
}
if (length(within_orfs)==0){ # if not keep all of them
  gr_orfs_subset=gr_orfs
}

# define ensemble of potentially overlapping orfs (1 single, 2 or more)
ensembles = GenomicRanges::reduce(gr_orfs_subset, ignore.strand=TRUE, min.gapwidth=0)
n_ensembles=length(ensembles)
ensembles_outcome=list()

# for each ensemble of overlapping orfs (i index), resolve (minimize overlap and maximize coding density)

threshold_overlap=as.numeric(args[3]) # only orfs with less than this value of overlapp are authorized

for (i in 1:length(ensembles)){ # i index : an ensemble of n overlapping orfs
  print(paste0("ensemble name: ",ensembles[i]))
  print(paste0("   ensemble of overlapping orfs number: ",i, " out of ", n_ensembles))

  interval = ensembles[i]
  # identify corresponding orfs
  orfs_id_within_this_interval=queryHits(findOverlaps(gr_orfs_subset, interval, type = "any", ignore.strand=TRUE))
  orfs_within_this_interval=gr_orfs_subset[orfs_id_within_this_interval]

  # how many orfs ?
  norfs_within_this_ensembl=length(orfs_within_this_interval)
  print(paste0("   number of orfs within this ensemble : ", norfs_within_this_ensembl))

  ##############################################################
  # simple case : only a single orf in this "ensembl" of orfs  #
  ##############################################################
  if(norfs_within_this_ensembl==1){
    ensembles_outcome[[i]]=orfs_within_this_interval
  }

  ##############################################################
  # more complicated case : severals orfs in this ensembl :
  #   resolve by maximizing coding length  #
  ##############################################################

  if(norfs_within_this_ensembl>1){

    # evaluate all combination of orfs from 2 to norfs_within_this_ensembl
    combn_list=list()
    for (k in 1:(norfs_within_this_ensembl)){
      combn_list[[k]]=combn(x = c(1:norfs_within_this_ensembl), m = k)
    }

    # determine if the orfs are compatible or not (overlap or not) for each combination (=column)
    longest_compatible_orfs_combinations=lapply(X = combn_list, FUN=function(x){
      n_orf_sampled=dim(x)[1]
      n_combinations=dim(x)[2]
      res=as.data.frame(matrix(data = NA, nrow = n_orf_sampled+1,ncol = n_combinations))
      for (c in 1:n_combinations){ # for each combination test whether orfs overlap.
        # if not calculate orf length
        orfs=orfs_within_this_interval[x[,c]]
        res[1:n_orf_sampled,c]=overlapsAny(query = orfs, drop.self=T, minoverlap = threshold_overlap, ignore.strand=TRUE) # test whether they overlap
        if (sum(res[1:n_orf_sampled,c])==0) { # if =0 they are not overlaping
          total_orf_length=sum(width(orfs)) # calculate orfs length
          res[n_orf_sampled+1,c]=total_orf_length
        }
      }
      # choose best
      if (sum(is.na(res[n_orf_sampled+1,]))<n_combinations){
        best=which(res[n_orf_sampled+1,]==max(res[n_orf_sampled+1,], na.rm = T))
        return(x[,best])
      }
      else {
        return(NA)
      }
    })

    # determine the coding length in the different scenario (with 1, 2, 3 etc... orfs)
    orf_length_combinations=unlist(lapply(X = longest_compatible_orfs_combinations, FUN = function(x){
      if (is.na(x[1])) {
        return(NA)
      }
      else {
        orf_length=sum(width(orfs_within_this_interval[as.vector(x)]))
        return(orf_length)
      }
    }))

    winner_comb=which(orf_length_combinations==max(orf_length_combinations, na.rm=TRUE))

    ensembles_outcome[[i]]=orfs_within_this_interval[longest_compatible_orfs_combinations[[winner_comb]]]
  }
}


################################################@
# END  ORF PREDICTION FUNCTION
################################################@


final_table=lapply(X = ensembles_outcome, FUN=function(x){
  as.data.frame(x)
})
final_table=do.call(rbind.data.frame, final_table)
head(final_table)

# reshape as a gff file
source="getorf_JV"
feature="gene"
score="."
phase="."
attibutes=final_table$width

final_table$source=source
final_table$feature=feature
final_table$score=score
final_table$phase=phase
final_table$attibutes=attibutes
gff=final_table[, c(1,6, 7, 2,3,8,5, 9, 10)]
head(gff)

# add sequence length (usefull for plots)
t=data.frame(names(genome), width(genome))
gff2=merge(gff, t, by.x = "seqnames", by.y="names.genome.", all.x=TRUE, all.y=TRUE)

### Write gff file

outfile_gff=paste0(tools::file_path_sans_ext(args[1]), paste0("_prediction_option", args[4],".gff"))
write.table(gff, file = outfile_gff, row.names = F, col.names = F, quote = FALSE, sep = "\t")

# write the protein sequences
feat_seq=list()
gr_outcomes=NULL
gr_outcomes=ensembles_outcome[[1]]
for (i in 2:length(ensembles_outcome)){
  gr_outcomes=c(gr_outcomes,ensembles_outcome[[i]])
}

feat_seq <- BSgenome::getSeq(genome, gr_outcomes)
feat_seq <- translate(x = feat_seq)
gr_outcomesdf=as.data.frame(gr_outcomes)
names=NULL
for (i in 1:dim(gr_outcomesdf)[1]){
  names[i]=paste(gr_outcomesdf[i,1], gr_outcomesdf[i,2], gr_outcomesdf[i,3], gr_outcomesdf[i,5], sep="_")
}
names(feat_seq)=names

writeXStringSet(feat_seq, filepath = paste0(tools::file_path_sans_ext(args[1]), paste0("_prediction_option", args[4],".fa")))


# plot predictions
tab=gff2
levels(tab$strand)=c("TRUE", "FALSE", "NA")
contigs=unique(tab$seqnames)
n_contigs=length(contigs)
plots=list()
for (i in 1:n_contigs){
tab_subset=tab[which(tab$seqnames==contigs[i]),]
plots[[i]]=ggplot(tab_subset, aes(xmin = start, xmax = end, y = seqnames, forward = strand)) +
  geom_gene_arrow() + facet_wrap(~ seqnames, scales = "free_y", ncol = 1)  +
  theme_genes() + geom_segment(aes(y = seqnames, yend = seqnames, x=width.genome.), xend = 100000, colour = "white", size = 2) +
  xlim(0,tab_subset$width.genome.[1]) +xlab("")+ylab("")
}

ggsave(
  filename = paste0(tools::file_path_sans_ext(args[1]), "_prediction_option", args[4],".pdf"),
  plot = marrangeGrob(plots, nrow=1, ncol=1),
)
