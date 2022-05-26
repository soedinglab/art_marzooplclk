#R script to pull out only orthologs to sequences of interest (SOIs).
#These will then be fed to InterProScan for annotation.
#
#This is done by parsing the .tsv files in the Orthologues/
#subdirectory in OrthoFinder's output directory.
#
#At this stage, I am simply going to find the orthogroup IDs
#that correspond to SOIs, and put them in a table with the
#sample names, and the orthogroup identifiers.
#
#Then I will use this information to pull the necessary sequences
#out of inputs directory I created for the OrthoFinder run.
#(This is located in under newrun/outputs/orthofinder/inputs.)
#What I will do is parse the TSVs I have now to get a two
#column dataframe containing the file name and the transcript(s)
#name, and then write this out to one file per file name containing
#all the transcripts names from that file name.
#I will then feed these "lists" to seqkit grep (I have it installed
#via conda) to grep the inputs/ I mentioned earlier for the specific
#sequences.
#I am specifically using that inputs/ directory because it has all
#the FASTA files--including those of the references--in one spot.


#--------------------------------------------------------------------------------

rm(list = ls())


library(dplyr)
library(stringr)
library(magrittr)
library(tidyr)
library(data.table)


#--------------------------------------------------------------------------------

#Function that takes a directory path as input
#and returns a dataframe of all OrthoFinder ortholog TSV files
#as output, with an additional column for the file name.
ofdirstodf <- function(x){
  dat <- bind_rows(lapply(list.files(x, full.names = TRUE), function(x){
    dat <- fread(x)
    names(dat) <- c("orthogroup", "sp1", "sp2")
    dat$filename <- basename(x)
    return(dat)
  }))
  return(dat)
}


#--------------------------------------------------------------------------------


#Paths.
mypath <- "/path/to/outputs" #SET THIS PATH HERE.

#Identifying the directories containing the references
#and orthologs to sequences from the references.
refdirs <- list.files(mypath, pattern = "_UP.*_allrefs", full.names = TRUE, include.dirs = TRUE)


#Reading in the OrthoFinder data using the custom
#function from above.
ofdat <- bind_rows(lapply(refdirs, ofdirstodf))


#Filtering down to retain SOI matches only.
ofdat %<>% filter(str_detect(sp1, "SOIREF_"))

#Separating each transcript matched into its own row.
ofdat %<>% separate_rows(sp2, sep = ",")
ofdat %<>% mutate(sp2 = str_replace_all(sp2, "[\\s\\t]+", ""))

#Also moving the corresponding references into separate rows, as
#there can be many-to-many pairwise orthologs.
ofdat %<>% separate_rows(sp1, sep = ",")
ofdat %<>% mutate(sp1 = str_replace_all(sp1, "[\\s\\t]+", ""))

#Creating reference and target file names from the filename column.
ofdat %<>% separate(filename, into = c("refname", "tarname"), sep = "__v__")
#Mutating refname and tarname columns to prepare for output files.
ofdat %<>% mutate(tarname = str_replace(tarname, "\\.tsv", "_ofcands.txt"))
ofdat %<>% mutate(refname = paste0(refname, "_ofrefs.txt"))

#Unique sequences from the references.
soirefs <- ofdat %>% 
  select(c(refname, sp1)) %>% 
  group_by(refname) %>% 
  distinct(sp1, .keep_all = TRUE) %>%
  ungroup()

#Unique sequences from the candidates.
soicands <- ofdat %>% 
  select(c(tarname, sp2)) %>% 
  group_by(tarname) %>% 
  distinct(sp2, .keep_all = TRUE) %>%
  ungroup()





#There are two sets of sequences to be written.
#soirefs are the reference sequences that matches candidates. As many-to-many
#orthologs do exist in the dataset, this also contains some OTHREF_ sequences.
#soicands are all sequences (including OTHREFs and SOIREFs) matched to the soirefs.

#Path to OrthoFinder directory containing all the input files (including references).
sfpath <- "/cbscratch/vraghav/newrun/outputs/orthofinder/inputs"

#First grepping and writing out the soicands.

#Creating output directory if it doesn't exist already.
outpath <- "/cbscratch/vraghav/newrun/outputs/orthofinder/rscript_of_getsois_ofcands"
if(!dir.exists(outpath)){dir.create(outpath)}


for(i in unique(soicands$tarname)){
  
  #i <- unique(soicands$tarname)[[1]]
  
  cat("Writing out SOIs from ", i, "!!\n")
  curdat <- soicands %>% filter(tarname == i) %>% distinct(sp2)
  
  
  write.table(curdat, file = paste0(outpath, "/", i), sep = "\n", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  
  
  infilename <- paste0(sfpath, "/", str_replace(i, "_ofcands.txt", ".fasta"))
  soifilename <- paste0(outpath, "/", i)
  outfilename <- paste0(outpath, "/", str_replace(i, "_ofcands.txt", "_ofcands.fasta"))
  
  cat("Done!! Now grepping SOIs from ", basename(infilename), " and putting them in ", basename(outfilename), "!!\n")
  
  mycmd <- paste0("source /home/mpg08/vraghav/miniconda3/etc/profile.d/conda.sh; conda activate seqkit_conda; seqkit grep -f ", soifilename, " ", infilename, " -o ", outfilename, "; conda deactivate") 
  
  system(mycmd)
  
  cat("Done!!\n")
  
}



#Now grepping and writing out the soirefs.
#Putting the refs and cands in separate directories.
outpath <- "/cbscratch/vraghav/newrun/outputs/orthofinder/rscript_of_getsois_ofrefs"
if(!dir.exists(outpath)){dir.create(outpath)}

for(i in unique(soirefs$refname)){
  
  #i <- unique(soirefs$refname)[[1]]
  
  cat("Writing out SOIREFSs from ", i, "!!\n")
  curdat <- soirefs %>% filter(refname == i) %>% distinct(sp1)
  
  
  write.table(curdat, file = paste0(outpath, "/", i), sep = "\n", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  
  
  infilename <- paste0(sfpath, "/", str_replace(i, "_ofrefs.txt", ".fasta"))
  soifilename <- paste0(outpath, "/", i)
  outfilename <- paste0(outpath, "/", str_replace(i, "_ofrefs.txt", "_ofrefs.fasta"))
  
  cat("Done!! Now grepping SOIs from ", basename(infilename), " and putting them in ", basename(outfilename), "!!\n")
  
  mycmd <- paste0("source /home/mpg08/vraghav/miniconda3/etc/profile.d/conda.sh; conda activate seqkit_conda; seqkit grep -f ", soifilename, " ", infilename, " -o ", outfilename, "; conda deactivate") 
  
  system(mycmd)
  
  cat("Done!!\n")
  
}





#DONE.

#--------------------------------------------------------------------------------

