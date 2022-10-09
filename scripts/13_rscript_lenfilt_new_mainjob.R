#Creates a representative set of Trinity-TransDecoder proteins
#by selecting the longest protein for each Trinity "gene".


rm(list = ls())

mypath <- "/path/to/outputs" #SET THIS PATH HERE.


library(seqinr)
library(tidyr)
library(magrittr)
library(stringr)
library(dplyr)


#----

#Custom list_files() function that checks whether the listed "files" are directories, and 
#excludes directories automatically.
list_files <- function(path = ".", pattern = NULL, all.files = FALSE, full.names = TRUE, 
                       recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE, 
                       incl_dirs = FALSE){
  #path <- paste0(mypath, "/", refpath)
  if(path == ".") { path = getwd() }
  
  #Include directories if recursive is set.
  if(incl_dirs & recursive){include.dirs = TRUE}
  
  #Needs to have full.names = TRUE in order to get full path to pass to dir.exists()
  files <- list.files(path = path, pattern = pattern, all.files = all.files, full.names = TRUE, 
                      recursive = recursive, ignore.case = ignore.case, include.dirs = include.dirs, 
                      no.. = no..)
  
  if(!incl_dirs){
    files <- files[!dir.exists(files)]
  }
  
  if(!full.names){
    return(basename(files))
  } else{
    return(files)
  }
  
}
#Faster compiled version.
#list_files_cmp <- compiler::cmpfun(list_files)




#----



#Function to return a data.frame given the path to a fasta file.
#The data.frame contains the sequence's name, sequence, and the filename (in that order).
#This function needs the seqinr package.
fastodf <- function(x){
  
  require(seqinr)
  
  #refpath <- "rbh_ref_prep"
  #refpat <- "*.fasta"
  #x <- list.files(paste0(mypath, "/", refpath), pattern = refpat, full.names = TRUE)[5]
  curdat <- seqinr::read.fasta(file = x, as.string = TRUE, whole.header = TRUE, seqtype = "AA")
  
  #curdat
  #curseq <- unlist(curdat)
  #curname <- names(curseq)
  
  curdat <- data.frame(curname = names(unlist(curdat)), curseq = unlist(curdat), stringsAsFactors = FALSE)
  curdat$filename <- rep_len(basename(x), length.out = nrow(curdat))
  
  #Just to check if all the sequence names have been assigned properly.
  #all(row.names(curdat) == curdat$curname)
  
  row.names(curdat) <- NULL
  
  return(curdat)
  
}



#----


#Function to return a data.frame given the path to a directory containing a collection of fasta files.
#Calls fastodf internally. Needs the seqinr package.
fasdirdf <- function(path = NULL, pat = NULL){
  #tst <- 
  #path <- paste0(mypath, "/", refpath)
  #pat <- "*.fasta"
  if(is.null(path)){ path <- getwd() }
  if(is.null(pat)){ pat <- "*" }
  #return(do.call("rbind", lapply(list.files(path = path, pattern = pat, full.names = TRUE), fastodf)))
  return(do.call("rbind", lapply(list_files(path = path, pattern = pat, full.names = TRUE, incl_dirs = FALSE), fastodf)))
}




#----


#MAIN SCRIPT BODY. THE FUNCTIONS ABOVE ARE AUXILIARY FUNCTIONS NEEDED FOR THE
#MAIN BODY BELOW.



#inpath <- paste0(mypath, "/", "outputs", "/", "transdecoder_trinity")
inpath <- paste0(mypath, "/", "outputs", "/", "mmseqs2_protclu_td")
outpath <- paste0(mypath, "/", "outputs", "/", "rscript_lenfilt")

if(!dir.exists(outpath)){ dir.create(outpath, recursive = TRUE) }

#For testing
#inpath <- paste0(mypath, "/", "transdecoder_trinity")
#outpath <- paste0(mypath, "/", "rlenfilt")

#Creating list of input files.
#inpfiles <- list_files(path = inpath, pattern = "_fixed_noast.fasta")
inpfiles <- list_files(path = inpath, pattern = "_pclu_rep_seq.fasta")

#Looping through and creating the outputs.
for(i in 1:length(inpfiles)){
  
  #i <- 1
  curfile <- inpfiles[i]

  curdat <- fastodf(curfile)
  
  cat("Input file: ", curfile, " with ", nrow(curdat), " sequences including multiple isoforms!!\n")
  
  #Just in case there are some curname whitespace shenanigans.
  curdat %<>% mutate(curname = str_replace(curname, "[\\s\\t]+$", ""))
  
  #Calculating lengths, gene and ORF-type columns.
  curdat %<>%
    rowwise() %>%
    mutate(slen = nchar(curseq), 
           gene = str_extract(curname, "^.*(?=_i[0-9]+\\.)"), 
           orf_type = str_extract(curname, "(?<=_ORFtype_).*$"))
  curdat <- data.frame(curdat, stringsAsFactors = FALSE)
  
  #Creating ORF-complete/incomplete column. 1 if complete, -1 otherwise.
  curdat %<>% 
    mutate(comp_orf = ifelse(str_detect(orf_type, "complete"), 1, -1))
  
  #Selecting the longest isoform as the representative for each gene.
  #Note: there are, in many cases, multiple longest isoforms that are
  #like 90% identical. This is basically a toss up in this regard, since
  #all of them are equally valid.
  curdat %<>% 
    group_by(gene) %>%
    arrange(desc(slen), desc(comp_orf), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
  #  slice_max(slen, n = 1, with_ties = FALSE) %>%
  #  ungroup()
  
  
  #Preparing to write out the file.
  curdat %<>%
    transmute(seq = paste0(">", curname, "_len_", slen, "\n", curseq))
  
  #outfile <- str_replace(basename(curfile), "_transdecoder_fixed_noast.fasta$", "_td_rfilt.fasta")
  outfile <- str_replace(basename(curfile), "_pclu_rep_seq.fasta$", "_td_rfilt.fasta")
  
  write.table(curdat, file = paste0(outpath, "/", outfile), sep = "\n",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat("Output file: ", outfile, " with ", nrow(curdat), " sequences (one isoform per gene)!!\n")

  
}

cat("All done!!\n")


#----------------------------------------------------------------------------------------------------------------#

