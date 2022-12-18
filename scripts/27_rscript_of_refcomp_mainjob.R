#R script for comparing whether OrthoFinder's sensitivity is affected by the
#addition of extra reference proteomes in comparison to the run with just
#the three reference proteomes.


#--------------------------------------------------------------------------------


rm(list = ls())

#Base path for this working directory.
mypath <- "/path/to/of_refcomp" #SET THIS PATH HERE.
#All relevant data--and this script--are located here or in subdirectories.


#Setting working directory.
setwd(mypath)

#Sourcing custom functions used in this script.
#Path to this script needs to be set w.r.t. working directory.
#source("21a_rscript_ccsel_auxfunc.R")



#Libraries.
library(ggvenn) #For the Venn diagram.
library(data.table) #For fread().
library(purrr) #Required by custom function fast_fread().
library(magrittr)
library(tidyr)
library(stringr)
library(dplyr)



#--------------------------------------------------------------------------------


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



#--------------------------------------------------------------------------------


#Main outputs directory path.
mainoutdir <- paste0(mypath, "/", "outputs")


#Loading in proteins to look at data. This is a subset of all proteins I
#have tagged as being of interest.
selseqs_path <- paste0(mypath, "/", "cc_orthofinder_dbs", "/", "seqs_of_interest_final.csv")
selseqs <- data.table::fread(selseqs_path, header = TRUE)
rm(selseqs_path)



#Reading in the OrthoFinder orthologs data.
#Orthologues_ext <- data with extra reference species. (ofext)
#Orthologues_min <- data with only the 3 main reference species. (ofmin)

#Reading in ofext
ofpath <- paste0(mainoutdir, "/", "Orthologues_ext")
refdirs <- list.files(ofpath, pattern = "_UP.*_allrefs", 
                      full.names = TRUE, include.dirs = TRUE)
#Reading in the data using the custom ofdirstodf function from rscript_ccsel_auxfunc.R.
ofext <- bind_rows(lapply(refdirs, ofdirstodf))
#rm(ofpath, refdirs)

#Reading in ofmin
ofpath <- paste0(mainoutdir, "/", "Orthologues_min")
refdirs <- list.files(ofpath, pattern = "_UP.*_allrefs", 
                      full.names = TRUE, include.dirs = TRUE)
#Reading in the data using the custom ofdirstodf function from rscript_ccsel_auxfunc.R.
ofmin <- bind_rows(lapply(refdirs, ofdirstodf))
rm(ofpath, refdirs)


#Adding a column to indicate which comparison set the data is from
#before merging them into a single data.frame.
ofext %<>% mutate(src = "ofext")
ofmin %<>% mutate(src = "ofmin")

ofdat <- bind_rows(ofext, ofmin)
rm(ofext, ofmin)


#Filtering to retain only those rows where refseq contains a SOIREF.
#The column can contain more than one sequence at this point.
ofdat %<>% filter(str_detect(refseq, "SOIREF"))

#Putting each reference value and match value pairs in refseq into its own row.
ofdat %<>% separate_rows(refseq, sep = ",")
ofdat %<>% mutate(refseq = str_replace_all(refseq, "[\\s\\t]", ""))
ofdat %<>% separate_rows(matchseq, sep = ",")
ofdat %<>% mutate(matchseq = str_replace_all(matchseq, "[\\s\\t]", ""))

#I have a subset of all the SOIREFs I have that I want to focus on.
#I've loaded this into selseqs. Filtering down to this dataset.
selseqs %<>% select(c(short_id, upid)) %>% mutate(upid = str_remove(upid, "\\-.*$"))
ofdat %<>% mutate(upid = str_extract(refseq, "(?<=_)[A-Za-z0-9]+(?=_[a-z]+$)"))
ofdat %<>% filter(upid %in% selseqs$upid)


#Extracting species identifiers from the filenames.
ofdat %<>% mutate(spref = str_extract(filename, "^[A-Za-z]+_[a-z]+"), 
                  spmatch = str_extract(filename, "(?<=__v__)[A-Za-z]+_[a-z]+"))


#Filtering to retain only those rows that correspond
#to SOIREF in refseq and TRINITY in matchseq.
ofdat %<>% filter(str_detect(refseq, "^SOIREF") & str_detect(matchseq, "^TRINITY_"))

#Retaining only in-use columns.
ofdat %<>% select(-c(orthogroup, filename, upid, spref))

#All I want to do is see if the ofext and ofmins have the same matchseqs for their refseqs.
#So I will create a UID that consists of spmatch and refseq, and pivot wider on this.
ofdat %<>% mutate(uid = paste(spmatch, refseq, sep = "__"), .keep = "unused")


# #Splitting and fulljoining.
# ofext <- ofdat %>% filter(src == "ofext") %>% select(-src) %>% rename("matchseq_ofext" = matchseq)
# ofmin <- ofdat %>% filter(src == "ofmin") %>% select(-src) %>% rename("matchseq_ofmin" = matchseq)
# #Full joining to see which "new" matches (if any) have been discovered
# out <- full_join(ofext, ofmin, by = "uid") %>% select(c(uid, matchseq_ofmin, matchseq_ofext))
# out %>% summarize(n_min = sum(!is.na(matchseq_ofmin))/nrow(.), n_ext = sum(!is.na(matchseq_ofext))/nrow(.))


#For each matchseq, checking which OrthoFinder data set (min or ext, or both)
#it is present in.
ofdat %<>% mutate(pres = 1) %>% pivot_wider(names_from = src, values_from = pres, values_fill = 0)

#Checking how many are present in one but not the other.
#Assuming the union of both sets to be the ground truth, then
#the proportion that is missing in each can be calculated like so:
ofdat %>% summarize(n_min = sum(ofmin == 0), n_ext = sum(ofext == 0))

#Checking if these are predominantly complete or incomplete.
ofdat %<>% mutate(iscmp = ifelse(str_detect(matchseq, "ORFtype_complete"), TRUE, FALSE))
ofdat %>% summarize(n_min = sum(ofmin == 0 & iscmp), n_ext = sum(ofext == 0 & iscmp))


#Venn diagram.
ofmin <- ofdat %>% filter(ofmin == 1) %>% select(matchseq) %>% as.vector() %>% unlist()
names(ofmin) <- NULL
ofext <- ofdat %>% filter(ofext == 1) %>% select(matchseq) %>% as.vector() %>% unlist()
names(ofext) <- NULL
vennlist <- list("min" = ofmin, "ext" = ofext)
#Plot.
plt <- ofdat %>%
  mutate(across(starts_with("of"), as.logical)) %>%
  ggplot() +
  geom_venn(aes(A = ofext, B = ofmin),
            set_names = c("With additional reference proteomes (n = 137)", 
                          "No additional reference proteomes"), 
            text_size = 6, set_name_size = 6) + theme_void() #+ coord_fixed()
plt

outfile <- paste0(mainoutdir, "/", "of_refcomp_venndiag.pdf")
ggsave(filename = outfile, plot = plt, width = 20, height = 10, dpi = 600, device = cairo_pdf)

