#R script for identifying and extracting candidate orthologs to
#known circadian clock sequences.


#ENSURE THAT THE seqs_of_interest_final.csv FILE FROM THE
#scripts/of_ref_prep directory IS IN THE SAME DIRECTORY
#THIS SCRIPT IS IN.


#Inputs:
#OrthoFinder's Orthologues/sp1__v__sp2 files - only those files where
#the "query" (sp1) is a reference proteome.
#InterProScan - TSV files for all reference and candidate sequences as (identified
#by OrthoFinder)
#FASTA files of the reference and candidate sequences mentioned above.


#--------------------------------------------------------------------------------


rm(list = ls())
cat("R script to find orthologs using OrthoFinder outputs!!\n")



#Sourcing custom functions used in this script.
#Path to this script needs to be set w.r.t. working directory.
source("rscript_ccsel_auxfunc.R")

#Base path for this working directory.




mypath <- "/path/to/outputs" #SET THIS PATH HERE.



#All relevant data--and this script--are located here or in subdirectories.

#Setting working directory.
setwd(mypath)



#Libraries.
#library(Biostrings)
#library(protr)
library(seqinr) #For read.fasta() used by my custom FASTA reader.
#library(ape) #For the read.gff() function.
library(data.table) #For fread().
library(purrr) #Required by custom function fast_fread().
library(magrittr)
library(tidyr)
library(stringr)
library(dplyr)



#--------------------------------------------------------------------------------


#Main outputs directory path.
mainoutdir <- paste0(mypath, "/", "outputs")

#Basic approach:
#OrthoFinder identifies pairwise orthologs.
#The idea is to simply use this feature to extract all pairwise orthologs of
#the reference circadian clock sequences.
#Then to "confirm" the data, InterProScan data (e.g., PANTHER family affiliation)
#can be used.



#Loading the OrthoFinder pairwise orthologs data.
#Only thing that will need to be updated here upon a new OrthoFinder run is
#the Results_DATE directory path.
ofpath <- paste0(mainoutdir, "/", "orthofinder/outputs/Results_Dec13/Orthologues")
cat("Reading in the OrthoFinder pairwise matches data from ", ofpath, "!!\n")
#Identifying the directories containing the  "reference" sp1__v__sp2 TSV files.
refdirs <- list.files(ofpath, pattern = "_UP.*_allrefs", 
                      full.names = TRUE, include.dirs = TRUE)
#Reading in the data using the custom ofdirstodf function from rscript_ccsel_auxfunc.R.
ofdat <- bind_rows(lapply(refdirs, ofdirstodf))
rm(ofpath, refdirs)


#Reading in the InterProScan annotations.
#Ensure that the quote symbol passed to fast_fread is "". It seems to cause R
#to crash otherwise.
#REFERENCE SEQUENCES.
ipspath_refs <- paste0(mainoutdir, "/interproscan_ofrefs")
cat("Reading in InterProScan data for the OrthoFinder references from ", ipspath_refs, "!!\n")
ipsnames_refs <- c("prot_acc", "seq_md5_dig", "seq_len", "analysis_plat", 
                   "signature_acc", "signature_desc", "start_loc", "stop_loc", 
                   "match_score", "match_status", "date_of_run", "ipannot_acc", 
                   "ipannot_desc", "go_annot", "pathway_annot")
ipsrefs <- fast_fread(locpath = ipspath_refs, 
                      locpattern = "_res_fixed.tsv$", locnamelist = ipsnames_refs, 
                      mysep = "\t", myquote = "")
rm(ipspath_refs, ipsnames_refs)
#CANDIDATE ORTHOLOG SEQUENCES.
ipspath_cands <- paste0(mainoutdir, "/interproscan_ofcands")
cat("Reading in InterProScan data for the OrthoFinder candidates from ", ipspath_cands, "!!\n")
ipsnames_cands <- c("prot_acc", "seq_md5_dig", "seq_len", "analysis_plat", 
                    "signature_acc", "signature_desc", "start_loc", "stop_loc", 
                    "match_score", "match_status", "date_of_run", "ipannot_acc", 
                    "ipannot_desc", "go_annot", "pathway_annot")
ipscands <- fast_fread(locpath = ipspath_cands, 
                       locpattern = "_res_fixed.tsv$", 
                       locnamelist = ipsnames_cands, mysep = "\t", myquote = "")
rm(ipspath_cands, ipsnames_cands)



#Reading in the FASTA files.
#All reference sequences + orthogroup OTHREF_ sequences on the query side from OrthoFinder.
faspath_refs <- paste0(mainoutdir, "/orthofinder/rscript_of_getsois_ofrefs")
cat("Reading in FASTA sequences for the references from ", faspath_refs, "!!\n")
fasrefs <- fasdirdf(path = faspath_refs, pat = ".fasta")
rm(faspath_refs)

#All matches to these references (+ orthogroup OTHREFs) above on the match side from OrthoFinder.
faspath_cands <- paste0(mainoutdir, "/orthofinder/rscript_of_getsois_ofcands")
cat("Reading in FASTA sequences for the candidates from ", faspath_cands, "!!\n")
fascands <- fasdirdf(path = faspath_cands, pat = ".fasta")
rm(faspath_cands)


#Although I've tentatively identified a large number of reference proteins, I am only 
#going to look at a select subset of those.
#I am loading the file of those identifiers below.
#These are located in the cc_orthofinder_dbs/ directory that should be at the same
#level as the outputs/ directory.
selseqs_path <- paste0(mypath, "/", "cc_orthofinder_dbs", "/", "seqs_of_interest_final.csv")
cat("Reading in the TSV file with the subset of the true sequences of interest from", selseqs_path, "!!\n")
selseqs <- data.table::fread(selseqs_path, header = TRUE)
rm(selseqs_path)
#selseqs %<>% select(c())


#Loading MMseqs2 annotation data.
cat("Loading MMseqs2 annotation data!!\n")
#This is from the filtered annotation data over at rscript_annotscmp/annotations.
mmseqs_path <- paste0(mainoutdir, "/", "rscript_annotscmp/annotations")
mmseqs_names <- c("curname", "sampname", "transcript", "trans_len", "prot_name", 
                  "prot_len", "orf_type", "swissprot_id", "swissprot_subjprotname", 
                  "swissprot_subjsp", "swissprot_gene", "swissprot_percid", "swissprot_evalue", 
                  "swissp_bitscore", "swissprot_qlen", "swissprot_slen", "swissprot_alnlen", 
                  "swissprot_qcov", "swissprot_tcov", "swissprot_subject", "eggnog_seed_ortholog", 
                  "eggnog_evalue", "eggnog_score", "eggNOG_OGs", "eggnog_max_annot_lvl", 
                  "eggnog_COG_category", "eggnog_Description", "eggnog_Preferred_name", "eggnog_GOs", 
                  "eggnog_EC", "eggnog_KEGG_ko", "eggnog_KEGG_Pathway", "eggnog_KEGG_Module", 
                  "eggnog_KEGG_Reaction", "eggnog_KEGG_rclass", "eggnog_BRITE", "eggnog_KEGG_TC", "eggnog_CAZy", 
                  "eggnog_BiGG_Reaction", "eggnog_PFAMs")
mmseqsdat <- fast_fread(locpath = mmseqs_path, 
                        locpattern = "_annotations.tsv$", 
                        locnamelist = mmseqs_names, mysep = "\t")
#Only need a couple of these columns.
mmseqsdat %<>% 
  select(c(filename, curname, swissprot_subject, 
           swissprot_percid, swissprot_evalue))
#mmseqsdat %<>% 
#  group_by(filename, curname) %>%
#  arrange(swissprot_evalue, desc(swissprot_percid), desc(swissprot_tcov), desc(swissprot_qcov), .by_group = TRUE) %>%
#  #slice_min(., order_by = swissprot_evalue, n = 1, with_ties = TRUE) %>%
#  slice_head(n = 1) %>%
#  ungroup()
#mmseqsdat %<>% 
#  select(c(curname, swissprot_subject, 
#           swissprot_percid, swissprot_evalue, 
#           filename))

rm(mmseqs_path, mmseqs_names)

#--------------------------------------------------------------------------------


#MAIN DATA PROCESSING SECTION


cat("Working with OrthoFinder data and filtering down to selected references and their candidates!!\n")

#Filtering ofdat to retain matches to reference SOIs only.
ofdat %<>% filter(str_detect(refseq, "SOIREF_"))
#Only doing this with the refseq, and not with the matchseq column because the
#pairwise orthologs are commutative.

#Looking to make sure there are no shenanigans in the data, references-wise.
test <- ofdat %>% distinct(refseq)
test %<>% mutate(line = row_number())
test %<>% separate_rows(refseq, sep = ", ")
test %<>% group_by(refseq) %>% mutate(lineset = paste0(line, collapse = " ; ")) %>% ungroup()
test %<>% distinct(refseq)
#Distinct brings the count down to as many sequence as fasrefs has. 
#The "duplicates" are because in some cases, pairwise orthologs had been defined
#at a deeper internal node for some match or the other, while a pairwise ortholog at a 
#lower (i.e., closer to the leaves) level also existed that included one of the same
#references.
#So everything is good here.
rm(test)



#Both refseq and matchseq columns have multiple items per row in some cases.
#This is because in some cases, the ortholog relationship is simply many-to-many.
#For instance TRINITY_DN8584_c0_g1_i1.p1_ORFtype_complete_len_713 from
#Phialella quadrata has a one-to-many orthology relationship with
#OTHREF_sp|P53762|ARNT_MOUSE, OTHREF_sp|Q61324|ARNT2_MOUSE, 
#SOIREF_ARNTL-2_Mus_musculus_Q2VPD4_main, SOIREF_ARNTL-1_Mus_musculus_Q9WTL8_main.
#This is simply because all these Mus paralogs evolved after the speciation event
#relating Mus and Phialella.
#I think I am happy accepting a number of false positives at this point.



#Both refseq and matchseq columns have multiple items in them in some cases.
#Moving these into their own rows. (This redundancy will be reduced later.)
ofdat %<>% separate_rows(matchseq, sep = ",")
ofdat %<>% mutate(matchseq = str_replace_all(matchseq, "[\\s\\t]", ""))
ofdat %<>% separate_rows(refseq, sep = ",")
ofdat %<>% mutate(refseq = str_replace_all(refseq, "[\\s\\t]", ""))


#Categorizing the refseq-matchseq combinations.
#SOIREF + TRINITY.
#SOIREF + SOIREF.
#SOIREF + OTHREF.
#OTHREF + TRINITY.
#OTHREF + SOIREF.
#OTHREF + OTHREF.
#Assigning by simply extracting the first \\w+ from the
#refseq and matchseq columns respectively.
ofdat %<>% mutate(refmatcat = paste0(str_extract(refseq, "^[A-Z]+(?=_)"), "__", 
                                     str_extract(matchseq, "^[A-Z]+(?=_)")))
#Not all the categories listed above would exist, as I am not importing the full
#symmetric set of pairwise orthologs here (i.e., if sp1__v__sp2 is here as a file,
#sp2__v__sp1 isn't as it makes no sense to just have all the matches listed twice
#over).
ofdat %>% distinct(refmatcat)


#Creating a reference species and match species column.
ofdat %<>% separate(filename, into = c("spref", "spmatch"), sep = "__v__")
ofdat %<>% mutate(spref = str_extract(spref, "^[A-Za-z]+_[A-Za-z0-0]+"), 
                  spmatch = str_extract(spmatch, "^[A-Za-z]+_[A-Za-z0-0]+"))


#--------------------------------



#For confirming the orthology of the matches, and filitering, all I need
#are the PANTHER assignments from InterProScan. PANTHER assigns the sequences
#to families (~ orthogroups) and subfamilies (functionally equivalent subsets of
#the orthogroups). So I will get the PANTHER data ready first.

#Since I need to do this extraction for both the reference and candidate sequences
#I wrapped the code into a function and placed it into rscript_ccsel_auxfunc.R.

cat("Working with PANTHER data from InterProScan, using for filtering!!\n")

#PANTHER data for the candidate sequences.
pancands <- ipspanext(ipscands)

#Preparing to merge this data into ofdat.
pancands %<>% 
  mutate(mergecol = paste0(filename, "__", prot_acc))
pancands %<>% select(c(mergecol, signature_acc, signature_desc))
names(pancands) <- c("mergecol", "pacc_match", "pdesc_match")

#Merging
ofdat %<>% mutate(mergecol = paste0(spmatch, "__", matchseq))
ofdat <- base::merge(ofdat, pancands, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)
ofdat <- data.frame(ofdat, stringsAsFactors = FALSE)
rm(pancands)



#Doing the same for the references.
panrefs <- ipspanext(ipsrefs)

#Preparing to merge this data into ofdat.
panrefs %<>% 
  mutate(mergecol = paste0(filename, "__", prot_acc))
panrefs %<>% select(c(mergecol, signature_acc, signature_desc))
names(panrefs) <- c("mergecol", "pacc_ref", "pdesc_ref")

#Merging
ofdat %<>% mutate(mergecol = paste0(spref, "__", refseq))
ofdat <- base::merge(ofdat, panrefs, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)
ofdat <- data.frame(ofdat, stringsAsFactors = FALSE)
rm(panrefs)

#Are any reference PANTHER annotations missing somehow?
ofdat %>% filter(is.na(pacc_ref))
#No.



#--------------------------------

#I will integrate the pfam data in a similar manner.
#I have a function written for this rscript_ccsel_auxfunc.R that can handle this
#for Pfam, CDD, and some other InterProScan databases.

cat("Working with Pfam data from InterProScan (not using for filtering)!!\n")

#First for the candidate sequences.
pfcands <- ipsdbext(ipscands, db = "Pfam")

#Preparing for merging
pfcands %<>% mutate(mergecol = paste0(filename, "__", prot_acc))
#Not retaining the match_score since I am not filtering by evalue.
pfcands %<>% select(c(mergecol, signature_acc, signature_desc, ssloc))

names(pfcands) <- c("mergecol", "pfacc_match", "pfdesc_match", "pfloc_match")

#Merging
ofdat %<>% mutate(mergecol = paste0(spmatch, "__", matchseq))
ofdat <- base::merge(ofdat, pfcands, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)
ofdat <- data.frame(ofdat, stringsAsFactors = FALSE)
rm(pfcands)



#Doing the same with the references.
pfrefs <- ipsdbext(ipsrefs, db = "Pfam")

#Preparing for merging
pfrefs %<>% mutate(mergecol = paste0(filename, "__", prot_acc))
#Not retaining the match_score since I am not filtering by evalue.
pfrefs %<>% select(c(mergecol, signature_acc, signature_desc, ssloc))

names(pfrefs) <- c("mergecol", "pfacc_ref", "pfdesc_ref", "pfloc_ref")

#Merging
ofdat %<>% mutate(mergecol = paste0(spref, "__", refseq))
ofdat <- base::merge(ofdat, pfrefs, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)
ofdat <- data.frame(ofdat, stringsAsFactors = FALSE)
rm(pfrefs)

#
rm(ipscands, ipsrefs)


#--------------------------------

bkup <- ofdat

#Candidate filtering.

cat("Performing candidate filtering!!\n")

#Some observations from manual data sifting.
#It looks like the Pfam domains sets CAN match totally, but the PANTHER assignments
#can be entirely different families. For example, TRINITY_DN7841_c0_g1_i1.p1_ORFtype_complete_len_302
#from Magelona mirabilis was picked up as an ortholog of mouse NR1D1 and NR1D2 (these are in-paralogs).
#Although the Pfam assignments for this match and its are identical in terms of domain compositions
#The sequences have been assigned to completely different PANTHER families.
#I think PANTHER assignments and groupings are fairly accurate, which means that these are just two
#different families that happen to share the same multi-domain architecture w.r.t. Pfam.

#SO I AM CHOOSING TO USE PANTHER AS MY CO-CLASSIFIER.

#Before everything else.
#Creating a protein category column.
ofdat %<>% mutate(protcat = str_extract(refseq, "(?<=^SOIREF_)[A-Za-z0-9\\-]+"))


#I only want to keep matches that have BOTH Pfam AND PANTHER annotations.
#No PANTHER means no classification confirmation.
#No Pfam means no domains = IMHO spurious false alignment.
ofdat %<>% filter(!is.na(pfacc_match) & !is.na(pacc_match))


#I will not retain matches that do not at least belong to the same PANTHER family.
ofdat %<>% 
  mutate(panfamref = str_extract(pacc_ref, "PTHR[0-9]+"), 
         panfammatch = str_extract(pacc_match, "PTHR[0-9]+")) %>%
  filter(panfamref == panfammatch) %>%
  select(-c(panfamref, panfammatch))


#One more thing I want to do is this.
#For each reference, there can be more than one match from the same species.
#For example, PDP1e has 6 matches from Calanus_helgolandicus.
#In this specific case, for example, I can see that there are actually two sequences
#that belong to the exact subfamily as the reference, and the others do not.
#Here I can actually filter out these "extra" matches, as it's unlikely that they are
#actually functionally equivalent to PDP1e (although they are orthologs).

#I will implement this by adding a "confidence" column. If the match and the reference
#belong to the same subfamily the confidence column will be 1, else 0.
ofdat %<>% mutate(mconf = ifelse(pacc_ref == pacc_match, 1, 0))

#Now to filter, all I need to do is group by the species and the reference sequence,
#and just keep all rows that have the maximum value for mconf. This will automatically
#keep only the exact PANTHER subfamily if the exist, and retain non-specific matches
#otherwise.
ofdat %<>% 
  group_by(spmatch, refseq) %>% 
  filter(mconf == max(mconf)) %>% 
  ungroup()


#--------------------------------


#I have a subset of all the SOIREFs I have that I want to focus on.
#I've loaded this into selseqs. Filtering down to this dataset.
selseqs %<>% select(c(short_id, upid)) %>% mutate(upid = str_remove(upid, "\\-.*$"))
ofdat %<>% mutate(upid = str_extract(refseq, "(?<=_)[A-Za-z0-9]+(?=_[a-z]+$)"))
ofdat %<>% filter(upid %in% selseqs$upid)
#rm(selseqs)
#ofdat %<>% select(-upid)


#Checking if all SOIs are in there (can be that something is missing simply
#because no matches were found)
foundprots <- ofdat %>% select(c(upid, protcat)) %>% distinct(upid, .keep_all = TRUE)
misprots <- selseqs %>% filter(!(upid %in% foundprots$upid))
#I will need this misprots table later to fill in for the final results table where
#I present which sequence have been found in which samples.


#Performing an initial check to see how many candidates per clock protein I have
#in each sample.
ofdat %>% 
  group_by(spmatch) %>% 
  count(protcat) %>% 
  pivot_wider(names_from = protcat, values_from = n)
#Okay, this looks quite acceptable.
#I do have situations where there is more than one candidate--e.g., with Crangon_crangon
#CLK, but this is to be expected. For instance, in that particular case there is a 3'
#sequence and an internal ORF. I COULD merge them, but there is no guarantee that those
#two fragments originated from the same parent gene. So I will retain the candidates as they
#are.

#Some additional checks.

#Have any matches occurred more than once?
ofdat %>% 
  group_by(spmatch) %>%
  count(matchseq) %>%
  filter(n > 1)
#No.

#Has any match been found by more than one reference?
ofdat %>%
  group_by(spmatch, matchseq) %>%
  count(protcat) %>%
  ungroup() %>%
  filter(n > 1)
#No.

rm(misprots, foundprots)


#Looking at how many matches have been found per protein category per species.

ofdat %>% 
  group_by(spmatch) %>% 
  count(protcat) %>% 
  ungroup() %>% 
  pivot_wider(names_from = protcat, values_from = n)

#This looks acceptable.
#Stopping filtering here, and moving on to output-related preparation.




#--------------------------------




cat("Incorporating the MMseqs2 functional annotation data (not for filtering)!!\n")

#Incorporating the mmseqs functional annotation data as a column.
#Preparing mmseqsdat for merging.
mmseqsdat %<>% 
  mutate(filename = str_replace(basename(filename), "_annotations.tsv$", ""))
mmseqsdat %<>% 
  mutate(mergecol = paste0(filename, "__", curname)) %>%
  select(-c(filename, curname))

#For this mergecol I need to remove the "grp" suffixes that exist in
#some of the sample names, as the annotations data have this sanitized.
ofdat %<>% 
  mutate(mergecol = paste0(str_replace(spmatch, "grp", ""), "__", matchseq))

#Filtering mmseqsdat to retain only sequences of interest before merging
#(Doing this just because those files are so huge.)
mmseqsdat %<>% 
  filter(mergecol %in% ofdat$mergecol)

#Merging.
ofdat <- merge(ofdat, mmseqsdat, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)
rm(mmseqsdat)



#--------------------------------




cat("Incorporating the FASTA sequences (not for filtering)!!\n")

#I will first incorporate the reference and match sequences themselves into the
#ofdat data.frame.

#Candidate sequences.
fascands %<>% mutate(filename = str_extract(filename, "^[A-Za-z]+_[a-z0-9]+"))
#I don't know why there's a space after the sequence names here, but there is.
fascands %<>% mutate(curname = str_replace_all(curname, "[\\s\\t]+", ""))
#For the reference side "candidate" matches in here (i.e., pairwise orthologs
#from the other references to the reference in question), I need to munge the
#sequence names as they contain the full sequence name string.
#For the SOIREFs in this set.
fascands %<>% mutate(curname = ifelse(str_detect(curname, "SOIREF_"), 
                                      str_replace_all(curname, "(?<=_main|_aux).*$", ""), curname))
#Also need to do these for the OTHREFs in this set.
fascands %<>% mutate(curname = ifelse(str_detect(curname, "OTHREF_"), 
                                      str_extract(curname, "^.*\\|[A-Z0-9]+_[A-Z]{5}"), curname))

fascands %<>% mutate(mergecol = paste0(filename, "__", curname))
fascands %<>% select(c(mergecol, curseq))
names(fascands) <- c("mergecol", "fas_match")

#Merging.
ofdat %<>% mutate(mergecol = paste0(spmatch, "__", matchseq))
ofdat <- base::merge(ofdat, fascands, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)


#Reference sequences.
fasrefs %<>% mutate(filename = str_extract(filename, "^[A-Za-z]+_[a-z0-9]+"))
#I don't know why there's a space after the sequence names here, but there is.
fasrefs %<>% mutate(curname = str_replace_all(curname, "[\\s\\t]+", ""))
#Munging the sequence names to match what's in ofdat.
fasrefs %<>% mutate(curname = ifelse(str_detect(curname, "SOIREF_"), 
                                     str_replace_all(curname, "(?<=_main|_aux).*$", ""), curname))
#Also need to do these for the OTHREFs in this set.
fasrefs %<>% mutate(curname = ifelse(str_detect(curname, "OTHREF_"), 
                                     str_extract(curname, "^.*\\|[A-Z0-9]+_[A-Z]{5}"), curname))

fasrefs %<>% mutate(mergecol = paste0(filename, "__", curname))
fasrefs %<>% select(c(mergecol, curseq))
names(fasrefs) <- c("mergecol", "fas_ref")

#Merging.
ofdat %<>% mutate(mergecol = paste0(spref, "__", refseq))
ofdat <- base::merge(ofdat, fasrefs, by = "mergecol", all.x = TRUE)
ofdat %<>% select(-mergecol)

rm(fasrefs, fascands)

#tst <- ofdat %>% filter(is.na(fas_match))


#--------------------------------

cat("Calculating pairwise similarity scores between the references and their matched candidates!!\n")

#I'd also like to add a column for the pairwise sequence similarity scores between
#the reference and the match.

#Sequence identity using Biostrings::pid() and Biostrings::pairwiseAlignment()
ofdat %<>% 
  rowwise() %>% 
  mutate(pid = Biostrings::pid(Biostrings::pairwiseAlignment(fas_ref, fas_match)))
ofdat <- data.frame(ofdat, stringsAsFactors = FALSE)

#Sequence similarity using protr::twoSeqSim(); this is quite fast.
#But the other option above should suffice, as more importantly, it provides me the
#percentage identity between the sequences.
#library(protr)
#ofdat %<>% 
#  rowwise() %>%
#  mutate(sim_score = protr::twoSeqSim(fas_ref, fas_match)@score)
#ofdat <- data.frame(ofdat, stringsAsFactors = FALSE)



#--------------------------------



#OUTPUTS SECTION.

#I will not do any visualizations, etc., within this script.
#So I want to write out data in such a way that it can be easily used for a broad
#variety of purposes.


#Setting up the output directory.
outdir <- paste0(mainoutdir, "/", "rscript_ccsel")
if(!dir.exists(outdir)){ dir.create(outdir) }

#Also need to fix the grp suffix in some species names.
ofdat %<>% mutate(spmatch = str_replace(spmatch, "grp$", ""))


#I need to add a column to ofdat that will carry the output name for the 
#candidate proteins, e.g., for the FASTA files.
#I'm going to make the name as follows:
#species + " " + clock protein name + " " + Sequence identifier.
ofdat %<>% mutate(candname = paste0(spmatch, "__", protcat, "__", matchseq))


#Arrange ofdat so that it is neat when written out.
ofdat %<>% arrange(protcat, spmatch)


#Will write two output tables with the full data.
#One for use with visualization programs and so forth.
#And the other for the publication supplement.

#Table for the publication supplement.
names(ofdat)
#Will not include mconf and upid columns.
#Columns not in use in the comments below. This is because
#the Pfam annotation columns have been taken out, as annotations
#will come in from elsewhere.
#pfloc_ref, pfloc_match,
ofdat_pub <- ofdat %>% 
  select(c(orthogroup, protcat, refmatcat, spref, spmatch, refseq, matchseq, 
           pid, pacc_ref, pdesc_ref, pacc_match, pdesc_match, 
           pfacc_ref, pfdesc_ref, pfacc_match, pfdesc_match,
           swissprot_subject, swissprot_evalue, swissprot_percid,
           fas_ref, fas_match, candname))
#Checking.
names(ofdat)[!(names(ofdat) %in% names(ofdat_pub))]


#Writing out the main table.
cat("Writing results to file!!\n")
#Rearranging the columns quickly.
maintab <- paste0(outdir, "/", "cc_cand_sel_main_table.csv")
cat("Writing main results table to ", maintab, "!!\n")
write.table(ofdat, maintab, sep = ",", 
            quote = TRUE, row.names = FALSE)

#Writing out publication table.
cat("Writing out publication table to file!!\n")
#Rearranging the columns quickly.
maintab <- paste0(outdir, "/", "cc_cand_sel_pub_table.csv")
cat("Writing main results table to ", maintab, "!!\n")
write.table(ofdat_pub, maintab, sep = ",", 
            quote = TRUE, row.names = FALSE)


#Writing out a per protein type FASTA file.
fasout <- paste0(outdir, "/", "fas_by_type")
cat("Writing out FASTA sequences per protein category to ", fasout, "!!\n")
if(!dir.exists(fasout)){ dir.create(fasout) }

for(i in 1:length(unique(ofdat$protcat))){
  
  #i <- 1
  
  curcat <- unique(ofdat$protcat)[i]
  outfile <- paste0("all_cands_", curcat, ".fasta")
  
  #Don't want to write out the references MATCHES into this as well.
  #So filtering those out here while selecting the right category.
  curdat <- ofdat %>% 
    filter(protcat == curcat) %>% 
    filter(!(spmatch %in% c("Drosophila_melanogaster", "Danaus_plexippus", "Mus_musculus"))) %>%
    select(c(refseq, fas_ref, candname, fas_match))
  
  #Matches can also have SOIREFs among them. To distinguish these from the actual reference
  #for the category, only the actual reference will have the FASTA header beginning with
  #SOIREF. All matches will begin the match's species name instead (as designed above in the
  #candname column).
  
  curdat %<>% 
    transmute(ref = paste0(">", refseq, "\n", fas_ref), 
              cand = paste0(">", candname, "\n", fas_match))
  
  curdat %<>% 
    pivot_longer(cols = everything(), names_to = "stype", values_to = "seq") %>%
    distinct(seq)
  
  #Writing out the table.
  write.table(curdat, file = paste0(fasout, "/", outfile), quote = FALSE, 
              row.names = FALSE, col.names = FALSE, sep = "\n")
}
rm(fasout, i, curcat, curdat, outfile)




#Writing out a per sample FASTA file.
fasout <- paste0(outdir, "/", "fas_by_samp")
cat("Writing out FASTA sequences per sample to ", fasout, "!!\n")
if(!dir.exists(fasout)){ dir.create(fasout) }

for(i in 1:length(unique(ofdat$spmatch))){
  
  #i <- 1
  
  curcat <- unique(ofdat$spmatch)[i]
  outfile <- paste0("all_cands_", curcat, ".fasta")
  
  curdat <- ofdat %>% 
    filter(spmatch == curcat) %>% 
    select(c(candname, fas_match))
  
  #Since the sequences are being grouped on the basis of samples here
  #no other filtering at this step is necessary.
  
  curdat %<>% 
    transmute(cand = paste0(">", candname, "\n", fas_match))
  
  #Writing out the table.
  #For the references, write out the files, but indicate that they are references.
  if((curcat %in% c("Drosophila_melanogaster", "Danaus_plexippus", "Mus_musculus"))){
    outfile <- paste0("REF_", outfile)
  }
  write.table(curdat, file = paste0(fasout, "/", outfile), quote = FALSE, 
              row.names = FALSE, col.names = FALSE, sep = "\n")
  
}
rm(fasout, i, curcat, curdat, outfile)


cat("All done!!\n")

#--------------------------------

#DONE.
