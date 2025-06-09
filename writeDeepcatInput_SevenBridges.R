# rmf 6.9.2025

#################
### FUNCTIONS ###
#################

writeDeepcatInputSB <- function(contigs, Tcells, outbase){
  
  # get list of barcodes corresponding to this sample
  print("Subsetting T cell object...")
  sub_Tcells <- subset(Tcells, subset = orig.ident == outbase)
  sub_Tcells_barcodes <- gsub("_.","",row.names(sub_Tcells[[]])) # remove underscore and suffix
  # ^ don't need this line if there aren't integrated samples in the object
  # but I don't think it will hurt anything to keep it here
  
  # keep only TRB locus rows for this sample and filter columns
  print("Filtering contigs...")
  contigs_filtered <- contigs %>%
    filter(cell_id %in% sub_Tcells_barcodes) %>%  # get T cells for just this sample
    filter(locus == "TRB") %>%  # filter out TCR alpha/delta/gamma/etc
    filter(junction != "" & junction_aa != "") %>%  # filter out cells without cdr3 nt or aa info
    select(consensus_count,  # number of reads for this contig
           cell_id, # cell barcode I think?
           junction,  # nucleotide sequence
           junction_aa,  # amino acid sequence
           v_call,  # V gene segment identified for this contig
           d_call,  # D gene segment identified for this contig
           j_call) %>% # J gene segment identified for this contig
    mutate(cdr3_length = nchar(junction))
  print(head(contigs_filtered))
  
  # calculate frequency of the clone
  print("Calculating frequency...")
  tmp <- contigs_filtered %>%
    select(cell_id, junction) %>%
    group_by(junction) %>%
    summarize(n_barcodes = n()) %>%
    arrange(desc(n_barcodes)) %>%
    mutate(freq = (n_barcodes/sum(n_barcodes))*100 ) %>%
    select(junction, freq)
  print(head(tmp))
  
  # merge the two dataframes to add frequency info
  print("Merging dataframes...")
  df <- merge(contigs_filtered, tmp, by = "junction", all.x = TRUE)
  
  # write to table
  write.table(contigs_filtered, file = paste("contig_info_",outbase,".txt",sep=""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # prepare for deepcat by filtering and renaming columns
  print("Prepping data for deepcat...")
  df <- df %>%
    rename(c("nucleotide" = "junction", "aminoAcid" = "junction_aa", "count (reads)" = "consensus_count",
             "frequencyCount (%)" = "freq", "cdr3Length" = "cdr3_length",
             "vMaxResolved" = "v_call", "dMaxResolved" = "d_call", "jMaxResolved" = "j_call"))
  print(head(df))
  
  ### add missing columns for DeepCAT input ###
  print("Adding v_gene columns...")
  newvec <- c()
  for (name in df$vMaxResolved){
    # remove everything after "*" in entries with "-" (indicates you can get the gene name)
    if (grepl("-",name) == FALSE){
      newname <- "unresolved"
    } else {
      newname <- gsub("\\*.*","",name)
    }
    newvec <- c(newvec, newname)
  }
  
  df$vGeneName <- newvec
  df$vFamilyName <- gsub("[-|\\*].*","",df$vMaxResolved)
  #df$vFamilyName <- gsub("-.*","",df$vMaxResolved)
  
  print("Adding d_gene columns...")
  newvec <- c()
  for (name in df$dMaxResolved){
    # remove everything after "*" in entries with "-" (indicates you can get the gene name)
    if (grepl("-",name) == FALSE){
      newname <- "unresolved"
    } else {
      newname <- gsub("\\*.*","",name)
    }
    newvec <- c(newvec, newname)
  }
  
  df$dGeneName <- newvec
  df$dFamilyName <- gsub("[-|\\*].*","",df$dMaxResolved)
  #df$dFamilyName <- gsub("-.*","",df$dMaxResolved)
  
  print("Adding j_gene columns...")
  newvec <- c()
  for (name in df$jMaxResolved){
    # remove everything after "*" in entries with "-" (indicates you can get the gene name)
    if (grepl("-",name) == FALSE){
      newname <- "unresolved"
    } else {
      newname <- gsub("\\*.*","",name)
    }
    newvec <- c(newvec, newname)
  }
  
  df$jGeneName <- newvec
  df$jFamilyName <- gsub("[-|\\*].*","",df$jMaxResolved)
  #df$jFamilyName <- gsub("-.*","",df$jMaxResolved)
  
  # required input for DeepCAT
  print("Adding extra columns for deepcat...")
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
  
  newcolnames <- setdiff(deepcat_colnames, colnames(df)) # get colnames not already in our dataframe
  
  df[newcolnames] <- NA  # add new colnames, fill with NA
  df <- df %>%
    select(c(cell_id, all_of(deepcat_colnames))) # reorder
  
  # write to table
  write.table(df, file = paste("deepcat_input_contigIDs_",outbase,".txt",sep=""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # write final input
  print("Writing input table...")
  df$cell_id <- NULL
  write.table(df, file = paste("deepcat_input_",outbase,".tsv",sep=""),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

############
### MAIN ###
############

library(Seurat)
library(tidyverse)

setwd("")

Tcells <- readRDS("")  # seurat object consisting of just T cells
contigs_file <- "" # VDJ_Dominant_Contigs_AIRR.tsv file
outbase <- "" # sample name corresponding to the orig.ident slot in the seurat object
# unique(Tcells$orig.ident)

contigs <- read.delim(contigs_file, sep = "\t")

# output is a tsv file with name format: paste("deepcat_input_",outbase,".tsv",sep="")
writeDeepcatInputSB(contigs, Tcells, outbase)
