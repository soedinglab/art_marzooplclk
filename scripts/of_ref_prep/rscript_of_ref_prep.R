#R script to prepare the reference proteomes for OrthoFinder.

#Reference proteomes are D. melanogaster, D. plexippus (just for CRY2), and M. musculus.
#All downloaded from UniProt.

#This script takes a CSV file indicating which sequences are of interest (SOIREFs), and manipulates
#the header strings to make them greppable later on. Also writes out sanitized FASTA files
#of the proteomes (and the SOIREFs) to the indicated output directories.

#MAKE SURE TO HAVE THE PROTEOMES IN A SUBDIRECTORY CALLED "proteomes",
#AND THE seqs_of_interest_reduced.csv file IN THIE SAME DIRECTORY THIS
#SCRIPT RESIDES IN.


#--------------------------------------------------------------------------------

rm(list = ls())

#Packages.

library("seqinr") #For reading in FASTA files
library("magrittr")
library("stringr")
library("dplyr")
library("tidyr")
library("data.table") #For fread()
library("httr") #For downloading external reference sequences

#Paths.
mypath <- "/path/to/ref_prep" #SET THIS PATH HERE. 



#--------------------------------------------------------------------------------

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

#Function to return a data.frame given the path to a fasta file.
#The data.frame contains the sequence's name, sequence, and the filename (in that order).
#This function needs the seqinr package.
fastodf <- function(x){
  
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



#Function to retrieve UniProt entries given a set of identifiers to a data.frame.
getsequp_todf <- function(upid = NULL, retfmt = "tab", retcols = c("id", "entry name", "protein names", "genes(PREFERRED)", "organism", "sequence", "reviewed")){
  
  if(is.null(upid)) { cat("No UniProt IDs provided!!\n"); break() }
  
  #retfmt = "tab"
  #retcols = c("id", "entry name", "genes(PREFERRED)", "organism", "sequence", "reviewed")
  
  baseurl <- "https://www.uniprot.org"
  resourceurl <- "uniprot"
  upurl <- paste0(baseurl, "/", resourceurl)
  
  retfmt <- paste0("&format=", retfmt)
  retcols <- paste0("&columns=", paste0(retcols, collapse = ","))
  retquer <- paste0("?query=accession:")
  
  #upid <- c("Q3TQ03", "E9GW67", "Q0QWP3")
  
  updat <- lapply(upid, function(x){
    
    #x <- "Q3TQ03"
    
    x <- paste0(upurl, "/", retquer, x, retfmt, retcols)
    
    #This is some weird stuff here. GET doesn't even have this listed. But this doesn't work if as = "text" is used in content() (where it is supposed to be).
    ret <- httr::GET(x, as = "text")
    httr::warn_for_status(ret)
    retdat <- httr::content(ret)
    retdat <- read.table(text = retdat, sep = "\t", header = TRUE)
    
  })
  #Closing connections to clear out lingering opened connections from httr::GET.
  closeAllConnections()
  
  updat <- do.call("rbind", updat)
}

#Function to retrieve UniProt entries given a set of identifiers to individual FASTA files saved to given path (getwd() by default).
getsequp_tofas <- function(upid = NULL, mypath = getwd(), fileprefix = "uniprot_"){
  
  if(is.null(upid)) { cat("No UniProt IDs provided!!\n"); break() }
  
  baseurl <- "https://www.uniprot.org"
  resourceurl <- "uniprot"
  upurl <- paste0(baseurl, "/", resourceurl)
  
  retfmt <- ".fasta"
  
  
  lapply(upid, function(x){
    
    download.file(paste0(upurl, "/", x, retfmt), destfile = paste0(mypath, "/", fileprefix, x, retfmt))
    
  })
  
  
  
}


#--------------------------------------------------------------------------------



#Data manipulation section.


#Reading in the proteomes.
#refpath <- paste0(mypath, "/", "of_ref_prep")
refprotpath <- paste0(mypath, "/", "proteomes")
refpat <- "*.fasta$"
refdat <- fasdirdf(path = refprotpath, pat = refpat)
rm(refprotpath, refpat)

#Importing file containing identifiers for sequences of interest.
soipath <- paste0(mypath, "/", "seqs_of_interest_reduced.csv")
soidat <- read.table(soipath, sep = ",", header = TRUE)
rm(soipath)
#Mutating upid column to suppress the trailing -[0-9] expression.
soidat %<>% mutate(upid = str_replace(upid, "\\-.*$", ""))


#Munging refdat to create a seqid column.
refdat %<>% mutate(seqid = str_extract(curname, "(?<=\\|).*(?=\\|)"))

#Creating a species identifier column
refdat %<>% mutate(srcsp = str_extract(filename, "^[A-Za-z]+_[a-z0-9]+(?=_)"))


#Identifying sequences of interest in refdat.
refdat %<>% mutate(soi = ifelse(seqid %in% soidat$upid, "SOI", NA))
#Checking if all SOIREFs have been ID'd properly.
refdat %>% 
  count(soi) %>% 
  filter(soi == "SOI") %>% 
  {ifelse(.$n == nrow(soidat), "All SOIREFs accounted for!!", "Not all SOIREFs accounted for!!")}
#All SOI are accounted for.


#Preparing a compressed version of soidat that will be merged with refdat.
#The short_id column is my idenitifier setup.
#The circ_role column is a string I am using to disambiguate core clock 
#proteins from auxiliary ones.
soismall <- soidat %>% select(c(upid, short_id, circ_role))
names(soismall) <- c("seqid", "short_id", "circ_role")
#Merging in the SOI data.
refdat <- merge(refdat, soismall, by = "seqid", all.x = TRUE)
rm(soismall)


#Constructing the sequence identifier string that will be inserted as the FIRST string
#in the FASTA header. The original header will follow this. As most tools truncate
#the sequence headers at the first whitespace, this sequence identifier will contain
#no internal spaces.
#This identifier will look like so: SOI_<short_id>_species_uniqueID_<main/aux> (uniqueID can be seqid, 
#i.e., UniProt or FlyBase accession string).
#For non-SOI sequences, it'll just be NOSOI_curname.
refdat %<>% mutate(headstr = ifelse(!is.na(soi), 
                                    paste0(">SOIREF_", short_id, "_", srcsp, "_", seqid, "_", circ_role, " ", curname),
                                    paste0(">OTHREF_", curname)))



#--------------------------------------------------------------------------------


#Prepping to write the reference transcriptomes (now with SOIs identified)
#to disk.
refdat %<>% select(c(headstr, curseq, filename))

#Output directory.
outpath <- paste0(mypath, "/", "final_ref_orthofinder")
if(!dir.exists(outpath)){dir.create(outpath)}

for(i in 1:length(unique(refdat$filename))){
  
  #i <- 1
  
  curdat <- refdat %>% filter(filename == unique(refdat$filename)[[i]])
  
  curdat %<>% transmute(seq = paste0(headstr, "\n", curseq))
  
  #All sequences set.
  outname <- str_replace(unique(refdat$filename)[[i]], "\\.fasta", "_allrefs.fasta")
  #outname <- str_replace(unique(refdat$filename)[[i]], "UP[0-9]+_[0-9]+\\.fasta", ".fasta")
  write.table(curdat, paste0(outpath, "/", outname), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  #SOIs only (from that particular dataset).
  outname <- str_replace(unique(refdat$filename)[[i]], "\\.fasta", "_soirefs.fasta")
  #outname <- str_replace(unique(refdat$filename)[[i]], "UP[0-9]+_[0-9]+\\.fasta", "_soirefs.fasta")
  curdat %<>% filter(str_detect(seq, "^>SOIREF_"))
  write.table(curdat, paste0(outpath, "/", outname), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}


#Writing out a separate dataset of ALL SOIS ONLY.
curdat <- refdat %>% filter(str_detect(headstr, "^>SOIREF_"))
curdat %<>% transmute(seq = paste0(headstr, "\n", curseq))
outname <- "orthofinder_sois_only.fasta"
write.table(curdat, paste0(outpath, "/", outname), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

#--------------------------------------------------------------------------------