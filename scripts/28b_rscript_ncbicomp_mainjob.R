#Script to process the hits from the NCBI transcriptome
#and genome assemblies that were downloaded and searched
#against (using MMseqs2) by ncbicomp_mainjob_noslurm.bash
#with the 10 selected reference circadian clock proteins as
#queries.

rm(list = ls())

library(kableExtra) #For LaTeX table
library(taxize) #For NCBI taxonomy
library(data.table)
library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)

mypath <- "/path/to/outputs/rscript_ncbicomp"

mainoutdir <- paste0(mypath, "/", "outputs")
ncbioutdir <- paste0(mainoutdir, "/", "ncbicomp")

hitdat <- list.files(ncbioutdir, pattern = "*.tab$", full.names = TRUE)
hitdat <- bind_rows(lapply(hitdat, function(x){
  df <- data.table::fread(x)
  df$fname <- x
  return(df)
}))
names(hitdat) <- c("qheader", "theader", "pident", "evalue", 
                   "bits", "qlen", "tlen", "alnlen", "qcov", 
                   "tcov", "fname")


hitdat %<>% 
  mutate(fname = str_replace(basename(fname), "\\.tab", "")) %>% 
  rowwise() %>%
  mutate(spmatch = unlist(str_split(fname, "_ncbi_"))[[1]], 
         samp = unlist(str_split(fname, "_ncbi_"))[[2]])

hitdat <- data.frame(hitdat, stringsAsFactors = FALSE)  

#First retaining only best hit for each target dataset
#and each query sequence against that target data set.
hitdat %<>% 
  group_by(fname, qheader) %>%
  #slice_min(order_by = evalue, with_ties = TRUE) %>%
  arrange(evalue, desc(pident), desc(tcov), desc(qcov), desc(alnlen), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

#Since this is a one-way search, will throw out everything that's in the
#twilight zone in terms of coverage.
hitdat %<>% filter(tcov >= 0.20 & qcov >= 0.20)

#Preparing to merge the two data.frames.
#Extracting protein catagory.
hitdat %<>% mutate(protcat = str_extract(qheader, "(?<=SOIREF_)[A-Za-z0-9\\-]+"))

#The data searched is a mixture of genomes and transcriptomes.
#Preferentially retain hits from genomes should matches from both
#genomes and transcriptomes be available for each protein category.
hitdat %<>% mutate(src_isgenome = ifelse(str_detect(theader, "TSA\\:"), 0, 1))
hitdat %<>% 
  group_by(spmatch, protcat) %>%
  arrange(desc(src_isgenome), evalue, desc(pident), desc(tcov), desc(qcov), desc(alnlen), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()
  
  

#To filter out only those hits that I need (matches for
#clock proteins not found by my own pipeline) I will intersect this data
#with the results tables from rscript_ccsel.R.
ccdat <- paste0(mainoutdir, "/", "rscript_ccsel", "/", "cc_cand_sel_pub_table.csv")
ccdat <- data.table::fread(ccdat)
ccdat %<>% 
  filter(refmatcat == "SOIREF__TRINITY") %>%
  select(c(protcat, spmatch, candname))
ccdat %<>%
  group_by(spmatch, protcat) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = protcat, values_from = n) %>%
  pivot_longer(cols = -spmatch, names_to = "protcat", values_to = "n") %>%
  filter(is.na(n))

ccdat %<>% select(-n)


#Merging
ccdat %<>% left_join(hitdat, by = c("spmatch", "protcat"))

#Writing this file out as is to disk.
ccout <- "ncbicomp_results.csv"
ccout <- paste0(ncbioutdir, "/", ccout)
write.table(ccdat, file = ccout, sep = ",", quote = TRUE, row.names = FALSE, col.names = TRUE)


#Tabulating the findings.
#Need to differentiate between protcats that are absent because
#no matches were found from protcats that are absent because
#no transcriptome was available for the search.

#For this I first need a list of the downloaded genomes/transcriptomes
#from the directory.
#This is already available.
assemlist <- read.table(file = paste0(ncbioutdir, "/", "ncbi_downloaded_genomes_proteomes.txt"))
names(assemlist) <- c("datname")
assemlist %<>% mutate(datname = str_replace(datname, "\\.fasta", ""))
assemlist %<>% 
  rowwise() %>%
  mutate(spmatch = unlist(str_split(datname, "_ncbi_"))[[1]])
assemlist <- data.frame(assemlist, stringsAsFactors = FALSE)
#Collecting all the assemblies for each sample into a single line.
assemlist %<>% 
  group_by(spmatch) %>%
  mutate(datname = paste0(datname, collapse = ",")) %>%
  distinct(spmatch, .keep_all = TRUE) %>%
  ungroup()


#Merging this with ccdat
ccdat %<>% full_join(assemlist, by = "spmatch")
#Retaining only those species with a protcat
#These are the ones that have something missing from
#the Orthofinder run.
ccdat %<>% filter(!is.na(protcat))

#Given that datname indicates whether or not data was available
#and qheader indicates whether or not a match was available,
#we can distinguish between hits missing because of a lack of data
#vs. because of missing homology.
ccdat %<>% mutate(matchstat = case_when(
  is.na(datname) ~ "ND",
  is.na(qheader) & !is.na(datname) ~ "NH",
  TRUE ~ "YES"
))
#Y - "Candidate detected", NH - "Homology absent", ND - "No data on NCBI"

#Ordering the proteins properly.
ccdat$protcat <- factor(ccdat$protcat, levels = c("CLK", "CYC", "TIM", "PER", 
                                                      "CRY1", "CRY2", "PDP1e", "VRI", 
                                                      "REV-ERBa", "RORa"))
ccdat %<>% arrange(protcat)

#Pivoting into table.
ccdat %<>% select(c(spmatch, protcat, matchstat)) %>%
  pivot_wider(names_from = protcat, values_from = matchstat)

ccdat


#----------------------------------------------------------------

##SPECIES GROUPING BY PHYLUM##

#Reading in table of species names and NCBI taxonomy IDs I prepared manually.
taxdat <- read.table(paste0(mypath, "/", "species_tree", "/", "species_data.csv"), sep = ",", header = TRUE)

#Filtering to retain non-reference species only.
taxdat %<>% filter(!str_detect(organism, "^Danaus|^Drosophila|^Mus"))
taxdat %<>% mutate(ncbi_taxid = as.numeric(str_extract(ncbi_taxid, "\\d+")))

#Getting the taxonomy data from NCBI.
ncbidat <- classification(taxdat$ncbi_taxid, db = "ncbi")

#Getting the phylogenetic tree for this.
#Will continue with the tree later.
ncbitree <- class2tree(ncbidat, check = TRUE)

#Row-binding the data.frames within this list.
ncbidat <- do.call("rbind", ncbidat)
#Extracting the NCBI IDs into a column.
ncbidat$ncbi_id <- rownames(ncbidat)
rownames(ncbidat) <- NULL
#The rank order is embedded in this column.
ncbidat %<>% separate(ncbi_id, into = c("ncbi_taxid", "rord"), sep = "\\.")

#Will keep the following taxonomic ranks.
#kingdom, phylum, class/subclass, order, family, genus, species.
#Animalia > Arthropoda > Copepoda > Calanoida > Temoridae > Temora > Temora longicornis
ncbidat %<>% 
  filter((rank %in% c("kingdom", "phylum", "class", "subclass", 
                      "order", "family", "genus", "species")))

#Dropping ID column, and pivoting wider.
ncbidat %<>%
  group_by(ncbi_taxid) %>%
  arrange(rord, .by_group = TRUE) %>%
  ungroup() %>%
  select(-c(id, rord)) %>%
  pivot_wider(names_from = rank, values_from = name)

#Merging this with my sample data.
taxdat %<>% select(-c(comments))
taxdat <- merge(taxdat, ncbidat, by = "ncbi_taxid", all.x = TRUE)
rm(ncbidat)

#Mutating phylum into a compact 2 letter code
taxdat %<>% mutate(phylum = str_extract(phylum, "^[A-Za-z]{2}"))
taxdat %<>% select(c(organism, phylum))

#Ordering taxdat to match tree.
taxdat$organism <- factor(taxdat$organism, 
                          levels = c("Acartia tonsa", "Acartia clausii", "Calanus helgolandicus", 
                                     "Centropages hamatus", "Temora longicornis", "Crangon crangon", 
                                     "Corystes sp.", "Hyperia sp.", "Podon leuckartii", "Evadne nordmanni", 
                                     "Poecilochaetus sp.", "Magelona mirabilis", "Phoronis muelleri", 
                                     "Oikopleura dioica", "Asterias rubens", "Rathkea octopunctata", 
                                     "Phialella quadrata"))

names(taxdat) <- c("disp_org", "disp_phy")
taxdat %<>% arrange(disp_org)
taxdat %<>% mutate(ord = row_number())

#----------------------------------------------------------------


#Setting up ccdat for LaTeX table.
ccdat %<>% mutate(disp_org = str_replace(spmatch, "_", " "), 
                 disp_org = ifelse(str_detect(disp_org, "\\ssp$"), paste0(disp_org, "."), disp_org), 
                 .keep = "unused")

#Merging in the taxonomy data.
ccdat %<>% full_join(taxdat, by = "disp_org")

#Sorting into the species tree order.
ccdat %<>% arrange(ord)

#Replacing NAs with *
ccdat %<>% mutate(across(.cols = -c(disp_org, disp_phy, ord), .fns = ~ifelse(is.na(.x), "-", .x)))
ccdat <- tibble::tibble(ccdat)

#Final table setup.
ccdat %>%
  select(matches("^[A-Z]{1}[A-Za-z0-9\\-]+", ignore.case = FALSE))


#Footnote:
genft <- "YES - candidate sequence found; NH - candidate not found because no homolog detected; ND - candidate not found because no data was available on NCBI; \"-\" - candidate discovered already by our workflow. Circadian clock protein abbreviations: CLK (Circadian locomotor output cycles kaput), CYC (Cycle), PER (period), TIM1 (Timeless) CRY1 (Cryptochrome 1), CRY2 (Cryptochrome 2), PDP1e (PAR Domain Protein 1 epsilon), VRI (Vrille), REV-ERBa (Nuclear Receptor Subfamily 1 Group D Member 1/NR1D1), RORa (RAR Related Orphan Receptor A; NR1F1/Nuclear receptor subfamily 1 group F member 1). Phyla: Ar - Arthropoda; An - Annelida; Ph - Phoronida; Ch - Chordata; Ec - Echinodermata; Cn - Cnidaria."
#Writing to object
conftab_out <- ccdat %>%
  select(c(disp_phy, disp_org), matches("^[A-Z]{1}[A-Za-z0-9\\-]+", ignore.case = FALSE)) %>%
  dplyr::rename("Organism" = disp_org, "Phylum" = disp_phy) %>%
  #arrange(Organism) %>%
  kbl(., format = "latex", booktabs = T, align = "clccccccccccc", linesep = "") %>%
  column_spec(2, italic = TRUE) %>%
  add_header_above(c(" " = 2, "Arthropod" = 8, "Mammalian" = 2)) %>%
  row_spec(row = 0, italic = FALSE, bold = TRUE) %>%
  footnote(general = genft, threeparttable = T) #%>%
#landscape()
#Here's the core latex table.
conftab_out

#Using save_kable with keep_tex = TRUE as a proxy to write the table out.
#Does not actually create the PDF since the GWDG R installation is missing
#dependencies.
#But the tex file does get written. This is all that's needed.
cat(conftab_out, file = paste0(ncbioutdir, "/", "ncbicomp_cands_table.tex"))
save_kable(conftab_out, 
           file = paste0(ncbioutdir, "/", "ncbicomp_cands_table.pdf"), 
           keep_tex = TRUE, self_contained = TRUE)

#----------------------------------------------------------------
