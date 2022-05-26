#Utility functions used in rscript_ccsel_main.R


#----

#Function that takes a directory path as input
#and returns a dataframe of all OrthoFinder ortholog TSV files
#as output, with an additional column for the file name.
ofdirstodf <- function(x){
  dat <- bind_rows(lapply(list.files(x, full.names = TRUE), function(x){
    dat <- fread(x)
    names(dat) <- c("orthogroup", "refseq", "matchseq")
    dat$filename <- basename(x)
    return(dat)
  }))
  return(dat)
}



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



#Fast file read function.
#For my version of default mmseqs column ordering (all columns?).
#But can be used with pretty much any data.
fast_fread <- function(locpath, locpattern = "*arth_ref_seqs_resultDB", locnamelist = NULL, mysep = "auto", myquote="\"", locheader = TRUE, ...){
  require("data.table")
  require("purrr")
  #locfiles <- list.files(path = locpath, pattern = locpattern, full.names = TRUE)
  rbindlist(lapply(list.files(path = locpath, pattern = locpattern, full.names = TRUE), function(x){
    locheader <- TRUE
    if(is.null(locnamelist)){
      locheader <- FALSE
    }
    if(!purrr::is_empty(data.table::fread(x, fill = T, sep = mysep, quote = myquote, header = locheader))){
      loctmpdf <- data.table::fread(x, stringsAsFactors = FALSE, fill = TRUE, sep = mysep, quote = myquote)
      if(isTRUE(locheader)){ names(loctmpdf) <- locnamelist }
      loctmpdf$filename <- x
    } else{
      if(isTRUE(locheader)){
        loctmpdf <- as.data.frame(matrix(ncol = length(locnamelist)))
        names(loctmpdf) <- locnamelist
        loctmpdf$filename <- x
      } else{
        #Do nothing
        cat(x, " has no data. Skipping!")
        #loctmpdf <- as.data.frame(matrix(NA))
        #loctmpdf$filename <- x
      }
    }
    return(loctmpdf)
  }))
}



#-----


#Compress and filter InterProScan data for a specific database. Cannot handle
#PANTHER (see function ipspanext instead).

ipsdbext <- function(x, db = NULL){
  #First for the candidate sequences.
  pf <- x %>% filter(analysis_plat == db)
  pf %<>% mutate(filename = str_extract(basename(filename), "^[A-Za-z]+_[a-z]+"))
  
  pf %<>% select(-c(seq_len, seq_md5_dig, match_status, date_of_run, 
                    ipannot_acc, ipannot_desc, go_annot, 
                    pathway_annot))
  pf %<>% mutate(signature_desc = str_replace_all(signature_desc, "[\\s\\t]+", "_"))
  
  pf %<>% 
    group_by(filename, prot_acc) %>% 
    arrange(start_loc, .by_group = TRUE) %>%
    ungroup()
  
  #Putting start and stop locations into a single value for each Pfam domain.
  pf %<>% 
    mutate(ssloc = paste0(start_loc, "__", stop_loc)) %>%
    select(-c(start_loc, stop_loc))
  
  #The domain calls are definitely duplicated, so getting that sorted out.
  pf %<>%
    group_by(filename, prot_acc) %>%
    distinct(signature_acc, .keep_all = TRUE) %>%
    ungroup()
  
  #Getting all the Pfam domains, start/stop locations, etc. per sequence into a single row.
  pf %<>%
    group_by(filename, prot_acc) %>%
    mutate(signature_acc = paste0(signature_acc, collapse = ";"),
           signature_desc = paste0(signature_desc, collapse = ";"), 
           match_score = paste0(match_score, collapse = ";"), 
           ssloc = paste0(ssloc, collapse = ";")) %>%
    distinct(signature_acc, .keep_all = TRUE) %>%
    ungroup()
  
  return(pf)
}


#----


#Function to extract and prepare PANTHER data given an InterProScan TSV table.
ipspanext <- function(x, db = "PANTHER"){
  
  pans <- x %>% filter(analysis_plat == db)
  pans %<>% mutate(filename = str_extract(basename(filename), "^[A-Za-z]+_[a-z]+"))
  
  pans %<>% select(-c(seq_len, seq_md5_dig, match_status, date_of_run, 
                      ipannot_acc, ipannot_desc, go_annot, 
                      pathway_annot, start_loc, stop_loc))
  pans %<>% mutate(signature_desc = str_replace_all(signature_desc, "[\\s\\t]+", "_"))
  
  
  #PANTHER annotation may be a family (PTHR1234) and/or a subfamily (PTHR1234:SF12).
  #I will assign these values to their own columns now.
  pans %<>% mutate(fam = str_extract(signature_acc, "PTHR[0-9]+"),
                   sfam = str_extract(signature_acc, "(?<=\\:)SF[0-9]+"))
  
  
  #I want to check to see if any sequences got assigned to more than one PANTHER family.
  #If THIS WERE TO BE THE CASE, I COULD ASSIGN EACH SEQUENCE TO ONE FAMILY ON 
  #THE BASIS OF E-VALUE.
  pans %>% 
    group_by(filename, prot_acc) %>%
    count(fam) %>%
    ungroup() #%>%
  #distinct(n)
  #No, so PANTHER assigns each sequence to a single family (or its sub-family below alongside)
  #But there are fully duplicated rows.
  
  #Removing fully duplicated rows for each unique sequence, if any.
  pans %<>% 
    group_by(filename, prot_acc) %>% 
    distinct(signature_acc, .keep_all = TRUE) %>%
    ungroup()
  
  #Now I will simply retain only the sub-family assignments if available, otherwise
  #the family identifier will be retained.
  pans %<>% 
    group_by(filename, prot_acc) %>%
    mutate(hassub = ifelse(str_detect(signature_acc, "\\:SF[0-9]+$"), 1, 0)) %>%
    arrange(desc(hassub)) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  return(pans) 
}


#----


