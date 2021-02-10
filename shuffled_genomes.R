
# library(tidyverse)
library(Biostrings)
library(seqinr)

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/genomes/"

setwd(dir)

fasta_files <- list.files(dir, pattern = ".fa", full.names = T)


shuffle_seqs <- function(fasta_file) {
  
  fs_out_name <- unlist(strsplit(basename(fasta_file),'[.]'))[1]
  
  x <- readDNAStringSet(fasta_file)
  
  # function to shuffle DNA sequences
  
  
  randomizeSeq <- function(mystring, kmer = 2){
    
    # myChars <- substring(mystring, 
    #           first = 1:(nchar(mystring) - 1), 
    #           last = kmer:nchar(mystring)-1)
    # split the string into characters
    myChars <- unlist(strsplit(as.character(mystring),''))
    # randomize the order of characters
    random.myChars <- sample(myChars,length(myChars))
    # concatenate them back into a string and return the result
    # DNAStringSet(paste(random.myChars,collapse=''))
    paste(random.myChars,collapse='')
  }
  
  # apply the randomize function to the DNAStringSet object 
  
  
  myRandomizedseqs <- sapply(x,randomizeSeq)
  myRandomizedseqs <- DNAStringSet(myRandomizedseqs)
  
  cat("save: ", fs_out_name, "\n at", getwd())
  
  writeXStringSet(myRandomizedseqs, paste0(fs_out_name, ".random.fs"))
  
  
  
}

# shuffle_seqs(fasta_files[1])
lapply(fasta_files, shuffle_seqs)
