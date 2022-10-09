#Script to go through the Trinity, rscript_lenfilt, diamond,
#and eggnogmapper files, and do the following things:

#Extract a "final transcriptome" dataset from the Trinity data
#using the rscript_lenfilt headers.

#Take all rscript_lenfilt names + diamond + eggnogmapper annotation
#files, and produce a single annotation file per sample.

rm(list = ls())

#Main path
mypath <- "/path/to/outputs" #SET THIS PATH HERE.

#Libraries

library("GO.db") #For annotated GO terms outputs table.
library(seqinr)
library(data.table)
library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)


#----------------------------------------------------------------
#Functions needed for this script.




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



#----


#Small function to calculate N50 for summary statistics.
n50_man <- function(x){
  ldf <- data.frame(lens = x)
  
  ldf %<>% 
    arrange(desc(lens)) %>% 
    mutate(csum = cumsum(lens)) %>%
    filter(csum >= round(max(csum)/2)) %>%
    filter(csum == min(csum)) %>%
    distinct(lens)
  
  return(as.numeric(ldf$lens))
}




#----------------------------------------------------------------


#This script will not read in the input data as a monolith at the
#beginning of the script as some of these datasets are really big.

#In fact, this script will run via a for loop.


#Outputs directory path.
mainoutdir <- paste0(mypath, "/", "outputs")

#I need a directory to put the final outputs in.
finaloutdir <- paste0(mainoutdir, "/", "rscript_annotscmp")
if(!dir.exists(finaloutdir)) { dir.create(finaloutdir) }

#Need output subdirectories.
#Transcriptomes.
outdir_trans <- paste0(finaloutdir, "/", "transcriptomes")
if(!dir.exists(outdir_trans)) { dir.create(outdir_trans) }
#Proteomes.
outdir_prots <- paste0(finaloutdir, "/", "proteomes")
if(!dir.exists(outdir_prots)) { dir.create(outdir_prots) }
#Annotations.
outdir_annots <- paste0(finaloutdir, "/", "annotations")
if(!dir.exists(outdir_annots)) { dir.create(outdir_annots) }
#eggnog-mapper orthologs.
outdir_emaporths <- paste0(finaloutdir, "/", "emapper_orthologs")
if(!dir.exists(outdir_emaporths)) { dir.create(outdir_emaporths) }
#GO terms table w/ annotations.
outdir_goterms <- paste0(finaloutdir, "/", "goterms")
if(!dir.exists(outdir_goterms)) { dir.create(outdir_goterms) }

#rscript_lenfilt's files form the crux of all operations here.
rlfiles <- list_files(paste0(mainoutdir, "/", "rscript_lenfilt"), pattern = "*.fasta")
rlfiles

#Path to eggnogmapper data.
emapdir <- paste0(mainoutdir, "/", "eggnogmapper")

#Path to diamond data.
#diadir <- paste0(mainoutdir, "/", "diamond")
diadir <- paste0(mainoutdir, "/", "mmseqs2_annots_rfilt")

#Path to trinity data.
trindir <- paste0(mainoutdir, "/", "trinity")


for(i in 1:length(rlfiles)){
  
  #i <- 1
  curfile <- rlfiles[[i]]
  
  #Need the sample name from the curfile string.
  cursname <- str_extract(basename(curfile), "[A-Za-z]+_[a-z]+")
  
  #If cursname has the grp identifier (indicating pooled libraries)
  #set up final name with that masked for outfiles.
  cursname_out <- str_replace(cursname, "grp$", "")
  
  cat("WORKING ON ANNOTATIONS AND FINAL SEQUENCE SETS FOR ", cursname_out, "!!\n")
  
  #I will create a directory for each sample under finaloutdir.
  #curoutdir <- paste0(finaloutdir, "/", cursname_out)
  #if(!dir.exists(curoutdir)) { dir.create(curoutdir) }
  
  #Reading in the rscript_lenfilt FASTA file.
  curlft <- fastodf(curfile)
  
  
  #----
  
  #SETTING UP THE ANNOTATION FILE.
  
  #Don't need the AA sequence, so dropping that column.
  curlft %<>% select(-curseq)
  
  #Quickly creating a transcript identifier column, sample name colulmn,
  #extracting ORF type and protein length columns, and renaming other columns.
  #Transcript ID
  curlft %<>% mutate(transcript = str_extract(curname, "^.*(?=\\.p[0-9]+)"))
  #Sample name
  curlft %<>% mutate(sampname = str_extract(filename, "[A-Za-z]+_[a-z]+"))
  #ORF type.
  curlft %<>% mutate(orf_type = str_extract(curname, "(?<=ORFtype_).*(?=_len)"))
  #Protein length.
  curlft %<>% mutate(prot_len = as.numeric(str_extract(curname, "(?<=_len_)[0-9]+$")))
  #Protein ID string.
  curlft %<>% mutate(prot_name = str_extract(curname, "^.*(?=_ORFtype)"))
  #Minimizing column set.
  curlft %<>% select(-filename)
  #curname column will stay on to act as key for merging eggnog and diamond
  #annotations.
  
  #----
  
  #READING IN TRANSCRIPTOME TO ADD DATA FROM IT TO ANNOTATION FILE.
  
  #Reading in the appropriate Trinity file.
  curtrin <- paste0(trindir, "/", cursname, "_trinity.fasta")
  
  cat("Incorporating annotatable data from Trinity transcriptome file ", curtrin, "!!\n")
  curtrin <- fastodf(curtrin)
  
  #I need the number of rows in curtrin as it indicates the number of transcripts in the
  #initial assembly. I need this for the summary statistics.
  #NUMBER OF TRANSCRIPTS IN INITAL ASSEMBLY.
  ntrans_raw <- nrow(curtrin)
  #Minimum transcript length, maximum, mean, N50.
  curtrin %<>% mutate(trans_len = nchar(curseq))
  trans_len_min_raw <- min(curtrin$trans_len)
  trans_len_max_raw <- max(curtrin$trans_len)
  trans_len_mean_raw <- mean(curtrin$trans_len)
  trans_n50_raw <- n50_man(curtrin$trans_len)
  #Dropping the trans_len column created above.
  curtrin %<>% select(-trans_len)
  
  #Creating a transcript ID column.
  curtrin %<>% mutate(transcript = str_extract(curname, "^[A-Za-z0-9\\_]+(?=\\s)"))
  
  #FILTERING TO RETAIN ONLY TRANSCRIPTS THAT ARE ALSO IN THE PROTEOME (I.E. FINAL ASSEMBLY)
  #Subsetting using this column.
  curtrin %<>% filter(transcript %in% curlft$transcript)
  
  #Calculating the transcript length. Never trust the values in the headers!!
  curtrin %<>% rowwise() %>% mutate(trans_len = nchar(curseq))
  curtrin <- data.frame(curtrin, stringsAsFactors = FALSE)
  
  #I just want to merge the transcript lengths with curlft.
  curtrin_mini <- curtrin %>% select(c(transcript, trans_len))
  curlft <- merge(curlft, curtrin_mini, by = "transcript", all.x = TRUE)
  rm(curtrin_mini)
  
  #Quickly arranging curlft columns.
  curlft %<>% select(c(curname, sampname, transcript, trans_len, 
                       prot_name, prot_len, orf_type))
  
  
  
  #----
  
  
  #WRITING OUT SANITIZED TRANSCRIPTOME AND PROTEOME FILES TO OUTPUT DIRECTORY.
  
  #Writing out the filtered transcriptome.
  cat("Munging Trinity FASTA file and writing it to directory ", outdir_trans, "!!\n")
  
  curtrin %<>% transmute(seq = paste0(">", transcript, "\n", curseq))
  curtrinout <- paste0(outdir_trans, "/", cursname_out, "_transcriptome.fasta")
  cat("Writing out the filtered transcriptome to ", curtrinout, "!!\n")
  write.table(curtrin, file = curtrinout, quote = FALSE, sep = "\n",
              row.names = FALSE, col.names = FALSE)
  rm(curtrin, curtrinout)
  
  cat("Writing transcriptome successful!!\n")
  
  #Taking this opportunity to also simply copy over curfile as cursname_proteome.fasta.
  cat("Also copying over proteome file to directory ", outdir_prots, "!!\n")
  protout <- paste0(outdir_prots, "/", cursname_out, "_proteome.fasta")
  cat("Copying over ", curfile, "as ", protout, "!!\n")
  #Will just use cp for this.
  syscmd <- paste0("cp ", curfile, " ", protout)
  execsys <- system(syscmd, intern = TRUE)
  if(length(execsys) == 0L) { cat("Copying proteome succesful!!\n")}
  rm(protout, syscmd, execsys)
  
  
  #----
  
  
  
  #CONTINUING SETTING UP ANNOTATIONS.
  
  
  
  #----
  
  #Need to incorporate MMseqs2 annotations into curlft.
  curmm <- paste0(diadir, "/", cursname, "_mmseqs2_outfmt6.tab")
  
  cat("Merging in MMseqs2 SwissProt annotations from ", curmm, "!!\n")
  
  curmm <- data.table::fread(curmm)
  names(curmm) <- c("curname", "swissprot_subject", "swissprot_percid", 
                    "swissprot_evalue", "swissp_bitscore", "swissprot_qlen", 
                    "swissprot_slen", "swissprot_alnlen", "swissprot_qcov", 
                    "swissprot_tcov")
  #First need to filter down to the best one match for each query by evalue and other params.
  curmm %<>%
    group_by(curname) %>%
    arrange(swissprot_evalue, desc(swissprot_percid), desc(swissprot_tcov), desc(swissprot_qcov), .by_group = TRUE) %>%
    #slice_min(., order_by = swissprot_evalue, n = 1, with_ties = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  #Extracting target species, target name, and target identifier.
  curmm %<>% mutate(swissprot_subjsp = str_extract(swissprot_subject, "(?<=OS\\=).*?(?=\\s[A-Z]+\\=)"))
  curmm %<>% 
    rowwise() %>%
    mutate(swissprot_subjprotname = str_match(swissprot_subject, 
                                              "(^sp\\|[A-Z0-9]+\\|[A-Z0-9_]+)\\s(.*(?=\\sOS\\=)).*$")[3], 
           swissprot_id = str_match(swissprot_subject, 
                                    "(^sp\\|[A-Z0-9]+\\|[A-Z0-9_]+)\\s(.*(?=\\sOS\\=)).*$")[2])
  curmm <- data.frame(curmm, stringsAsFactors = FALSE)
  #Can also extract the gene name.
  curmm %<>% mutate(swissprot_gene = str_extract(swissprot_subject, "(?<=GN\\=)[A-Za-z0-9\\-_]+"))
  
  #Rearranging columns.
  curmm %<>% select(c(curname, swissprot_id, swissprot_subjprotname, swissprot_subjsp, 
                      swissprot_gene, swissprot_percid, swissprot_evalue, swissp_bitscore, 
                      swissprot_qlen, swissprot_slen, swissprot_alnlen, swissprot_qcov,
                      swissprot_tcov, swissprot_subject))
  
  #Merging with curlft.
  curlft <- merge(curlft, curmm, by = "curname", all.x = TRUE)
  rm(curmm)
  
  cat("Done merging MMseqs2 SwissProt annotations!!\n")
  
  
  #----
  
  #SECTION NOT IN USE
  
  #INCORPORATING DIAMOND SWISSPROT ANNOTATIONS.
  
  
  #Need to incorporate Diamond annotations into curlft.
  #curdia <- paste0(diadir, "/", cursname, "_dmnd_outfmt6.tab")
  
  #cat("Merging in Diamond SwissProt annotations from ", curdia, "!!\n")
  
  #curdia <- data.table::fread(curdia, skip = 3)
  #names(curdia) <- c("curname", "swissprot_subject", "swissprot_percid", 
  #                   "swissprot_evalue", "swissp_bitscore", "swissprot_qlen", 
  #                   "swissprot_slen", "swissprot_alnlen")
  #Extracting target species, target name, and target identifier.
  #curdia %<>% mutate(swissprot_subjsp = str_extract(swissprot_subject, "(?<=OS\\=).*?(?=\\s[A-Z]+\\=)"))
  #curdia %<>% 
  #  rowwise() %>%
  #  mutate(swissprot_subjprotname = str_match(swissprot_subject, 
  #                                            "(^sp\\|[A-Z0-9]+\\|[A-Z0-9_]+)\\s(.*(?=\\sOS\\=)).*$")[3], 
  #         swissprot_id = str_match(swissprot_subject, 
  #                                  "(^sp\\|[A-Z0-9]+\\|[A-Z0-9_]+)\\s(.*(?=\\sOS\\=)).*$")[2])
  #curdia <- data.frame(curdia, stringsAsFactors = FALSE)
  #Can also extract the gene name.
  #curdia %<>% mutate(swissprot_gene = str_extract(swissprot_subject, "(?<=GN\\=)[A-Za-z0-9\\-_]+"))
  
  #Rearranging columns.
  #curdia %<>% select(c(curname, swissprot_id, swissprot_subjprotname, swissprot_subjsp, 
  #                     swissprot_gene, swissprot_percid, swissprot_evalue, swissp_bitscore, 
  #                     swissprot_qlen, swissprot_slen, swissprot_alnlen, swissprot_subject))
  
  #Merging with curlft.
  #curlft <- merge(curlft, curdia, by = "curname", all.x = TRUE)
  #rm(curdia)
  
  #cat("Done merging Diamond SwissProt annotations!!\n")
  
  #----
  
  
  #INCORPORATING EGGNOG-MAPPER ANNOTATIONS.
  
  #Reading in the eggnogmapper data.
  #eggnogmapper annotations.
  cureannots <- paste0(emapdir, "/", cursname, ".emapper.annotations")
  
  cat("Merging in eggnog-mapper annotations from ", cureannots, "!!\n")
  
  #Need to read this file in through the cmd argument and grepping out some columns
  #because it is not formatted in a standard way.
  cureannots <- data.table::fread(cmd = paste0("grep -v ", "\"^#\"", " ", cureannots), 
                                  sep = "\t", header = FALSE)
  names(cureannots) <- c("curname", "eggnog_seed_ortholog", "eggnog_evalue", "eggnog_score", "eggNOG_OGs", 
                         "eggnog_max_annot_lvl", "eggnog_COG_category", "eggnog_Description", 
                         "eggnog_Preferred_name", "eggnog_GOs", "eggnog_EC", "eggnog_KEGG_ko", 
                         "eggnog_KEGG_Pathway", "eggnog_KEGG_Module", "eggnog_KEGG_Reaction", 
                         "eggnog_KEGG_rclass", "eggnog_BRITE", "eggnog_KEGG_TC", "eggnog_CAZy", 
                         "eggnog_BiGG_Reaction", "eggnog_PFAMs")
  
  #Merging this with curlft
  curlft <- merge(curlft, cureannots, by = "curname", all.x = TRUE)
  rm(cureannots)
  
  cat("Done merging eggnog-mapper annotations!!\n")
  
  #----
  
  #COPYING OVER EGGNOGMAPPER ORTHOLOGS FILE TO OUTPUT DIRECTORY.
  
  #eggnogmapper orthologs.
  #This will not be merged into the table, but will simply be present in
  #the final directory.
  cat("Also copying over the eggnog-mapper orthologs file to ", outdir_emaporths, "!!\n")
  cureorth <- paste0(emapdir, "/", cursname, ".emapper.orthologs")
  eorthout <- paste0(outdir_emaporths, "/", cursname_out, "_emapper_orthologs.tsv")
  cat("Copying over ", cureorth, "as ", eorthout, "!!\n")
  syscmd <- paste0("cp ", cureorth, " ", eorthout)
  execsys <- system(syscmd, intern = TRUE)
  if(length(execsys) == 0L) { cat("Copying eggnog-mapper orthologs data succesful!!\n")}
  rm(eorthout, cureorth, syscmd, execsys)
  
  
  #----
  
  #Creating GO terms table.
  #Using functions from GO.db here.
  goterms <- curlft %>% 
    dplyr::select(c(sampname, curname, eggnog_GOs))
  goterms %<>% 
    separate_rows(eggnog_GOs, sep = ",")
  goterms %<>% filter(!is.na(eggnog_GOs))
  #Using functions from GO.db
  goterms %<>% 
    mutate(go_desc = Term(eggnog_GOs), go_cat = Ontology(eggnog_GOs),
           go_def = Definition(eggnog_GOs))
  #Filter out GO terms w/o go_desc.
  goterms %<>% filter(!is.na(go_desc))
  
  #Write this out to file.
  outfile <- paste0(outdir_goterms, "/", cursname_out, "_goterms.tsv")
  cat("Writing out the GO terms file to ", outfile, "!!\n")
  write.table(goterms, file = outfile, sep = "\t", quote = TRUE, row.names = FALSE)
  cat("Done writing GO terms annotations file!!\n")
  
  
  
  
  #----
  
  #WRITING OUT ANNOTATIONS FILE TO OUTPUT DIRECTORY.
  outfile <- paste0(outdir_annots, "/", cursname_out, "_annotations.tsv")
  cat("Writing out the main annotations file to ", outfile, "!!\n")
  write.table(curlft, file = outfile, sep = "\t", quote = TRUE, row.names = FALSE)
  cat("Done writing main annotations file!!\n")
  
  #----
  
  #Summary statistics for annotations.
  cat("Creating summary statistics for ", cursname_out, "!!\n")
  summtab_curr <- curlft %>% 
    summarise(sample = cursname_out, 
              nseqs_initial = ntrans_raw,
              n50_initial = trans_n50_raw,
              mean_translen_initial = trans_len_mean_raw,
              min_translen_initial = trans_len_min_raw,
              maxtranslen_initial = trans_len_max_raw,
              nseqs_final = nrow(.), 
              perc_cmp_final = sum(orf_type == "complete")/nrow(.),
              n50 = n50_man(trans_len),
              mean_translen = mean(trans_len, na.rm = TRUE), 
              min_translen = min(trans_len, na.rm = TRUE), 
              max_translen = max(trans_len, na.rm = TRUE),
              mean_protlen = mean(prot_len, na.rm = TRUE),
              min_protlen = min(prot_len, na.rm = TRUE), 
              max_protlen = max(prot_len, na.rm = TRUE), 
              annots_swissprot_only = sum(!is.na(swissprot_id) & is.na(eggnog_evalue)), 
              annots_eggnog_only = sum(is.na(swissprot_id) & !is.na(eggnog_evalue)), 
              annots_both = sum(!is.na(swissprot_id) & !is.na(eggnog_evalue)), 
              miss_annots = sum(is.na(swissprot_id) & is.na(eggnog_evalue)))
  
  if(i == 1){
    summtab <- summtab_curr
  } else{
    summtab <- bind_rows(summtab, summtab_curr)
  }
  
  
  #----
  
  
  rm(summtab_curr, outfile, curlft)
  cat("All operations done for ", cursname, "!!\n")
  cat("\n")
  
  
}

#Writing out the summary statistics table.
cat("Writing the summary statistics for all samples to ", finaloutdir, "!!\n")
summout <- paste0(finaloutdir, "/", "all_samps_", "summ_stats.tsv")
write.table(summtab, file = summout, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Written summary statistics to file ", summout, "!!\n")
rm(summtab, summout)

cat("Everything finished and done!!\n")


#Done.


#---------------------------------------------------------------------------------------------------------------------#

