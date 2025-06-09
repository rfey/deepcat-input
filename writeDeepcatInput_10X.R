# rmf 6.9.2025

#################
### FUNCTIONS ###
#################

writeDeepcatInput10X <- function(contigs, Tcells, clonotypes, outbase){
  # goal is to merge by clonotype id
  contigs_TRB <- contigs[contigs$chain == "TRB",] # filter out TCR alpha
  contigs_TRB <- contigs_TRB %>%
    select(barcode, length, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3_nt, reads, umis, raw_clonotype_id)
  
  # filter contigs: filter out barcodes that are not T cells
  Tcells <- subset(Tcells, subset = orig.ident == outbase)
  Tcell_barcodes <- gsub("_.","",row.names(Tcells[[]])) # remove underscore and suffix
  # ^ don't need this line if there aren't integrated samples in the object
  # but I don't think it will hurt anything to keep it here
  
  contigs_TRB <- contigs_TRB[contigs_TRB$barcode %in% Tcell_barcodes,]
  
  # filter clonotype info
  clonotypes_filtered <- clonotypes %>%
    select(c(clonotype_id, frequency, proportion))  # keep only necessary columns
  head(clonotypes_filtered)
  
  # merge metadata using clonotype id
  merged <- merge(contigs_TRB, 
                        clonotypes_filtered, 
                        by.x = "raw_clonotype_id",
                        by.y = "clonotype_id",
                        all.x = TRUE) # all.x=T keeps all extra rows in contigs
  head(merged)
  
  # write to table
  write.table(merged, file = paste("clonotype_barcodes_",outbase,".txt",sep=""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # prepare for deepcat
  merged <- merged %>%
    group_by(raw_clonotype_id) %>% 
    mutate(read_sum = sum(reads)) %>%
    distinct(cdr3, cdr3_nt, .keep_all = TRUE) %>% # removes rows with duplicate sequences, keeps all columns
    select(raw_clonotype_id, cdr3_nt, cdr3, read_sum, proportion, 
           length, v_gene, d_gene, j_gene) %>%  # filter and reorder columns
    rename(c("nucleotide" = "cdr3_nt", "aminoAcid" = "cdr3", "count (reads)" = "read_sum",
             "frequencyCount (%)" = "proportion", "cdr3Length" = "length", 
             "vMaxResolved" = "v_gene", "dMaxResolved" = "d_gene", "jMaxResolved" = "j_gene"))
  head(merged)
  
  ### add missing columns for DeepCAT input ###
  newvec <- c()
  for (name in merged$vMaxResolved){
    if (grepl("-",name) == FALSE){
      newname <- "unresolved"
    } else {
      newname <- gsub("\\*.*","",name)
    }
    newvec <- c(newvec, newname)
  }
  
  merged$vGeneName <- newvec
  merged$vFamilyName <- gsub("-.*","",merged$vMaxResolved)
  
  newvec <- c()
  for (name in merged$dMaxResolved){
    if (grepl("-",name) == FALSE){
      newname <- "unresolved"
    } else {
      newname <- gsub("\\*.*","",name)
    }
    newvec <- c(newvec, newname)
  }
  
  merged$dGeneName <- newvec
  merged$dFamilyName <- gsub("-.*","",merged$dMaxResolved)
  
  
  newvec <- c()
  for (name in merged$jMaxResolved){
    if (grepl("-",name) == FALSE){
      newname <- "unresolved"
    } else {
      newname <- gsub("\\*.*","",name)
    }
    newvec <- c(newvec, newname)
  }
  
  merged$jGeneName <- newvec
  merged$jFamilyName <- gsub("-.*","",merged$jMaxResolved)
  
  # required input for DeepCAT
  deepcat_colnames <- c("nucleotide","aminoAcid","count (reads)","frequencyCount (%)",
                        "cdr3Length","vMaxResolved","vFamilyName","vGeneName",
                        "vGeneAllele","vFamilyTies","vGeneNameTies", "vGeneAlleleTies", 
                        "dMaxResolved","dFamilyName","dGeneName","dGeneAllele",
                        "dFamilyTies", "dGeneNameTies","dGeneAlleleTies","jMaxResolved",
                        "jFamilyName","jGeneName","jGeneAllele", "jFamilyTies",
                        "jGeneNameTies","jGeneAlleleTies","vDeletion","n1Insertion",
                        "d5Deletion","d3Deletion","n2Insertion","jDeletion","vIndex",
                        "n1Index","dIndex","n2Index","jIndex","estimatedNumberGenomes",
                        "sequenceStatus","cloneResolved","vOrphon","dOrphon","jOrphon",
                        "vFunction","dFunction","jFunction","fractionNucleated")
  
  newcolnames <- setdiff(deepcat_colnames, colnames(merged)) 
  newcolnames
  
  merged[newcolnames] <- NA  # add new colnames, fill with NA
  merged <- merged %>%
    select(c(raw_clonotype_id, all_of(deepcat_colnames)))  # reorder
  
  # write to table
  write.table(merged, file = paste("clonotype_",outbase,".txt",sep=""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # write final input
  merged$raw_clonotype_id <- NULL
  write.table(merged, file = paste("deepcat_input_",outbase,".txt",sep=""),
              sep = "\t", quote = FALSE, row.names = FALSE)
}


############
### MAIN ###
############

library(Seurat)
library(tidyverse)

setwd("")

Tcells <- readRDS("")  # seurat object consisting of just T cells
contigs_file <- ""  # filtered_contig_annotations.csv
clonotypes_file <- ""  # clonotypes.csv
outbase <- "" # sample name corresponding to the orig.ident slot in the seurat object
# unique(Tcells$orig.ident)

contigs <- read.csv(contigs_file, check.names = FALSE)
clonotypes <- read.csv(clonotypes_file, check.names = FALSE)

# output is a tsv file with name format: paste("deepcat_input_",outbase,".tsv",sep="")
writeDeepcatInput(contigs, Tcells, clonotypes, outbase)

