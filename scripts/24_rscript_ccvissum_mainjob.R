#Script for visualizing and summarizing the candidate sequences found using
#rscript_ccsel_mainjob.R
#Must be run after that script.



rm(list = ls())

#Base path for this working directory.
#mypath <- "/home/mpg08/vraghav/rscript_ccvissum"
mypath <- "/home/owner/Nextcloud/laptop_rplace/eichele/art_marzooplclk"
#All relevant data--and this script--are located here.

#Setting working directory.
setwd(mypath)


#For phylogenetic ordering of samples during plotting
library(taxize)
#For ordering the sample names according to the phylogenetic
#tree order in the plots.
library(forcats)

#Uses the seqvisr package from GitHub, so install if necessary.
#devtools::install_github("vragh/seqvisr", force = TRUE)
library(ggtext) #For plotting
library(knitr) #For tables
library(kableExtra) #For tables
library(IRanges) #For protein domain selection from multiple overlapping domains
library(purrr)
library(seqvisr) #For the domain structure visualizations
library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(magrittr)
library(dplyr)

#--------------------------------------------------------------------------------


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



#--------------------------------------------------------------------------------




#

#Input data paths.
mainoutdir <- paste0(mypath, "/", "outputs")

#Reading in the main table from rscript_ccsel_main.R.
ccresdir <- paste0(mainoutdir, "/", "rscript_ccsel")
ccdf <- data.table::fread(paste0(ccresdir, "/", "cc_cand_sel_main_table.csv"))


#For constructing the domain structure diagrams, the
#InterProScan annotation tables will be used.
#Reading in the InterProScan annotations.
#Ensure that the quote symbol passed to fast_fread is "". It seems to cause R
#to crash otherwise.
#REFERENCE SEQUENCES.
ipspath_refs <- paste0(mainoutdir, "/interproscan_ofrefs")
ipsnames_refs <- c("prot_acc", "seq_md5_dig", "seq_len", "analysis_plat", 
                   "signature_acc", "signature_desc", "start_loc", "stop_loc", 
                   "match_score", "match_status", "date_of_run", "ipannot_acc", 
                   "ipannot_desc", "go_annot", "pathway_annot")
ipsrefs <- fast_fread(locpath = ipspath_refs, 
                      locpattern = "_res_fixed.tsv$", locnamelist = ipsnames_refs, 
                      mysep = "\t", myquote = "")
#Sites data.
ipsnames_refs_sites <- c("prot_acc", "seq_md5_dig", "seq_len", "analysis_plat", 
                   "signature_acc", "start_loc_dom", "stop_loc_dom", "n_res", "residue", 
                   "start_loc", "stop_loc", "signature_desc")
ipsrefs_sites <- fast_fread(locpath = ipspath_refs, 
                            locpattern = "_res.tsv.sites$", locnamelist = ipsnames_refs_sites, 
                            mysep = "\t", myquote = "")

rm(ipspath_refs, ipsnames_refs, ipsnames_refs_sites)

#CANDIDATE ORTHOLOG SEQUENCES.
ipspath_cands <- paste0(mainoutdir, "/interproscan_ofcands")
ipsnames_cands <- c("prot_acc", "seq_md5_dig", "seq_len", "analysis_plat", 
                    "signature_acc", "signature_desc", "start_loc", "stop_loc", 
                    "match_score", "match_status", "date_of_run", "ipannot_acc", 
                    "ipannot_desc", "go_annot", "pathway_annot")
ipscands <- fast_fread(locpath = ipspath_cands, 
                       locpattern = "_res_fixed.tsv$", 
                       locnamelist = ipsnames_cands, mysep = "\t", myquote = "")
#Sites data.
ipsnames_cands_sites <- c("prot_acc", "seq_md5_dig", "seq_len", "analysis_plat", 
                          "signature_acc", "start_loc_dom", "stop_loc_dom", "n_res", "residue", 
                          "start_loc", "stop_loc", "signature_desc")
ipscands_sites <- fast_fread(locpath = ipspath_cands, 
                            locpattern = "_res.tsv.sites$", locnamelist = ipsnames_cands_sites, 
                            mysep = "\t", myquote = "")

rm(ipspath_cands, ipsnames_cands, ipsnames_cands_sites)




#Also have othannots_ofrefs and othannots_ofcands directories
#containing sequence bias composition and NLS annotations.
#(fLPS and NLStradamus respectively.)
#oth_refs <- paste0(mypath, "/othannots_ofrefs")
#fLPS
#flps_refs <- bind_rows(lapply(list.files(oth_refs, pattern = "_flps.tab", full.names = TRUE), 
#                    function(x){df <- fread(x, skip = 4); df$filename <- x; return(df)}))
#names(flps_refs) <- c("seq", "bias", "lps", "start", 
#                "end", "res_count", "binom_p", "sig_desc", "filename")
#NLStradamus
#nls_refs <- fast_fread(locpath = oth_refs, locpattern =  "_nls.tab")
#names(nls_refs) <- c("seq", "algo", "score", "start", "stop", "actseq", "filename")
#rm(oth_refs)

#CANDIDATES.
#oth_cands <- paste0(mypath, "/othannots_ofcands")
#fLPS
#flps_cands <- bind_rows(lapply(list.files(oth_cands, pattern = "_flps.tab", full.names = TRUE), 
#                               function(x){df <- fread(x, skip = 4); df$filename <- x; return(df)}))
#names(flps_cands) <- c("seq", "bias", "lps", "start", 
#                       "end", "res_count", "binom_p", "sig_desc", "filename")
#NLStradamus
#nls_cands <- fast_fread(locpath = oth_cands, locpattern =  "_nls.tab")
#names(nls_cands) <- c("seq", "algo", "score", "start", "stop", "actseq", "filename")
#rm(oth_cands)

#Since these respective flps and nls annotations for candidates and references can be
#merged, doing so now.
#flps <- bind_rows(flps_refs, flps_cands)
#nls <- bind_rows(nls_refs, nls_cands)
#rm(flps_cands, flps_refs, nls_cands, nls_refs)


#Path to Pfam-A clans table.
#Need this to get the short identifiers for the Pfam domains for identification.
#pfamannot <- paste0(mypath, "/", "Pfam-A.clans.tsv")
#pfamannot <- data.table::fread(pfamannot, header = FALSE)



#--------------------------------------------------------------------------------

#Output directory for visualizations and stuff.
visoutdir <- "rscript_ccvissum"
visoutdir <- paste0(mainoutdir, "/", visoutdir)

if(!dir.exists(visoutdir)) { dir.create(visoutdir) }

#--------------------------------------------------------------------------------


##SPECIES GROUPING BY PHYLUM##

#Reading in table of species names and NCBI taxonomy IDs I prepared manually.
taxdat <- read.table(paste0(mainoutdir, "/", "species_tree", "/", "species_data.csv"), sep = ",", header = TRUE)

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


#Integrating phylogenetic orer data into visout
#visout %<>% mutate(disp_org = str_extract(seq, "^[A-Za-z]+\\s[a-z\\.]+"), 
#                  disp_org = ifelse(str_detect(disp_org, "\\ssp$"), paste0(disp_org, "."), disp_org)) %>%
#  full_join(taxdat, by = "disp_org") %>% 
#  mutate(seq = ifelse(str_detect(seq, "REF"), seq, paste0(seq, " (", disp_phy, ")")))
#Need to set the reference as ord = 0 to make sure it shows up on top.
#visout %<>% mutate(ord = ifelse(str_detect(seq, "REF"), 0, ord))

#----------------------------------------------------------------


#For confirming that the OrthoFinder approach actually works, I will use the
#reference vs. reference matches.

#Filtering from ccdf.
refconf <- ccdf %>% 
  filter(refmatcat == "SOIREF__SOIREF") %>%
  select(c(protcat, refseq, matchseq, spmatch))


#Notice NR1D1 is not here as it is has no ortholog in Drosophila or Danaus.
#So fill it in from a SOIREF__TRINITY match in ccdf.
#nr1d1 <- ccdf %>% filter(protcat == "NR1D1") %>% select(protcat, refseq, matchseq, spmatch) %>%
#  distinct(refseq, .keep_all = TRUE)
#nr1d1$matchseq <- NA
#nr1d1$spmatch <- NA
#Filling in.
#refconf <- bind_rows(refconf, nr1d1)
#rm(nr1d1)


#Extracting protein match category and UniProt accession.
refconf %<>% mutate(matchprotcat = str_extract(matchseq, "(?<=^SOIREF_)[A-Za-z0-9\\-]+"), 
                   matchupid = str_extract(matchseq, "[A-Z0-9]+(?=_main|_aux)"))
#And putting these into a single column.
refconf %<>% mutate(matchcat = paste0(matchprotcat, " ", "(", matchupid, ")"))
refconf %<>% select(-c(matchseq, matchprotcat, matchupid))

#Grouping the matchcat values for each protcat into single rows by species.
refconf %<>% 
  group_by(protcat, spmatch) %>% 
  arrange(matchcat, .by_group =  TRUE) %>% 
  mutate(matchcat = paste0(matchcat, collapse = ", ")) %>%
  ungroup() %>%
  distinct(matchcat, .keep_all = TRUE)

#Fixing match species names.
#I will use abbreviations throughout the table to make this easier.
#Fewer italics that way.
#refconf %<>% mutate(spmatch = str_replace(spmatch, "(?<=^[A-Z])[a-z]+_", ". "))
refconf$spmatch <- gsub("(^[A-Z]).*_([a-z]).*$", "\\1\\2", refconf$spmatch)

#Fixing that NA (NA) thing in NR1D1's matchcat cell.
#refconf$matchcat[refconf$protcat == "NR1D1"] <- NA

#Grouping each match species by protcat into a single row.
refconf %<>% 
  group_by(protcat, spmatch) %>%
  arrange(matchcat, .by_group = TRUE) %>%
  mutate(matchcat = paste0(matchcat, collapse = ", ")) %>%
  distinct(matchcat, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(matchcat = paste0(spmatch, ": ", matchcat)) %>%
  group_by(protcat) %>%
  arrange(spmatch, .by_group = TRUE) %>%
  mutate(matchcat = paste0(matchcat, collapse = "; ")) %>%
  distinct(matchcat, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-spmatch)

#I do want to have the matchcats by species in each row after all.
refconf %<>% separate_rows(matchcat, sep = "; ")

#Setting up the reference sequence column.
refconf %<>%
  mutate(refseq = str_replace_all(str_extract(refseq, "[A-Za-z]+_[a-z]+_[A-Z0-9]+"), "_", " "))
refconf$refseq <- gsub("(^[A-Z])[a-z]+\\s([a-z])[a-z]+\\s([A-Z0-9]+)$", 
                       "\\1\\2 \\3", 
                       refconf$refseq)

#Also adding in the protcat details here; protcat itself will stay on as a column
#to arrange the data, but will not be displayed in the table.
refconf %<>% 
  mutate(refseq = paste0(protcat, " ", "(", refseq, ")"))

#Adding a column to indicate which loop the proteins are a part of.
refconf %<>% mutate(ttfl_loop = case_when(
  str_detect(protcat, "CLK|CYC|TIM|PER|CRY2|CRY1|PDP1e|VRI") ~ "Arthropod",
  str_detect(protcat,  "REV\\-ERBa|RORa") ~ "Mammalian"
))

#Fixing that NA (NA) thing in NR1D1's matchcat cell.
#refconf$matchcat[refconf$protcat == "NR1D1"] <- NA


#Making protcat a factor column because I want the proteins grouped by
#the arm of the circadian cycle they belong to.
refconf$protcat <- factor(refconf$protcat, levels = c("CLK", "CYC", "TIM", "PER", 
                                                      "CRY1", "CRY2", "PDP1e", "VRI", 
                                                      "REV-ERBa", "RORa"))
refconf %<>% arrange(protcat)


#Rearranging the columns.
refconf %<>% select(c(protcat, ttfl_loop, refseq, matchcat))

#Quick glane.
refconf

#Building the latex table for export.
#I will handle the italicization manually.
genft <- "Species: Dm - Drosophila melanogaster, Dp - Danaus plexippus, Mm - Mus musculus. Protein name abbreviations: CLK (Circadian locomotor output cycles kaput), CYC (Cycle), PER (period), TIM (Timeless), CRY1 (Cryptochrome 1), CRY2 (Cryptochrome 2), PDP1e (PAR Domain Protein 1 epsilon), VRI (Vrille), REV-ERBa (NR1D1/Nuclear Receptor Subfamily 1 Group D Member 1), RORa (RAR Related Orphan Receptor A; NR1F1/Nuclear receptor subfamily 1 group F member 1), NPAS2 (Neuronal PAS Domain Protein 2), ARNTL (Aryl hydrocarbon Receptor Nuclear Translocator-Like protein), DBP (D-box Binding Protein), HLF (Hepatic Leukemia Factor), TEF (Thyrotroph Embryonic Factor), NFIL3 (Nuclear Factor, Interleukin 3 Regulated), EIP75B (Ecdysone-induced protein 75B, isoforms C/D; NR1D3), HR3 (Probable nuclear hormone receptor HR3). PER1, PER2, and PER3 are in-paralogs. ARNTL1 and ARNTL2 are in-paralogs."
refconf_tab <- refconf %>%
  select(-protcat) %>%
  dplyr::rename("TTFL" = ttfl_loop, "Reference" = refseq, "Ref. proteome match" = matchcat) %>%
  kbl(., format = "latex", booktabs = T, align = "lll", linesep = "") %>%
  #kable_styling(latex_options = c("scale_down"), position = "center", full_width = TRUE) %>%
  #kable_styling(position = "center", protect_latex = TRUE) %>%
  row_spec(row = 0, italic = FALSE, bold = TRUE) %>%
  collapse_rows(columns = 1:2, latex_hline = "none", valign = "top", row_group_label_position = "stack") %>%
  footnote(general = genft, threeparttable = T)

refconf_tab

#Using save_kable with keep_tex = TRUE as a proxy to write the table out.
#Does not actually create the PDF since the GWDG R installation is missing
#dependencies.
#But the tex file does get written. This is all that's needed.
cat(refconf_tab, file = paste0(visoutdir, "/", "of_ref_confirm.tex"))
save_kable(refconf_tab, 
           file = paste0(visoutdir, "/", "of_ref_confirm.pdf"), 
           keep_tex = TRUE, self_contained = TRUE)

rm(refconf, refconf_tab, genft)


#--------------------------------------------------------------------------------



#Visualizing the main matches table.
#Copying ccdf into a new data.frame.
confdf <- ccdf

#I can filter out all the non-SOIREF__TRINITY matches.
confdf %<>% filter(refmatcat == "SOIREF__TRINITY")

#Selecting columns.
#mconf is a column indicating whether the PANTHER annotation
#is identical or not.
#I am also calculating the ratios of the match sequence length
#to the reference sequence length here.
#And selecting mconf and this along with some other columns for
#creating the table.
confdf %<>% 
  rowwise() %>%
  mutate(reflen = nchar(fas_ref), matchlen = nchar(fas_match), 
         lenratio = as.numeric(matchlen)/as.numeric(reflen)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  select(c(protcat, spmatch, matchseq, mconf, pid, lenratio))


#Extracting ORFtype into a column.
confdf %<>%
  mutate(orftype = case_when(
    str_detect(matchseq, "_complete") ~ "cmp",
    str_detect(matchseq, "_internal") ~ "int",
    str_detect(matchseq, "_3prime") ~ "3p",
    str_detect(matchseq, "_5prime") ~ "5p"
  ))

confdf$orftype <- factor(confdf$orftype, levels = c("cmp", "3p", "5p", "int"))

#Counting the ORFtypes for each species and protein category combination.
#Also taking the opportunity to put the ORFtypes in the factor order from
#above.
confdf %<>% 
  group_by(protcat, spmatch) %>% 
  count(orftype) %>%
  arrange(orftype, .by_group = TRUE) %>%
  ungroup()

#Pivoting this wider.
confdf %<>% pivot_wider(names_from = orftype, values_from = n, values_fill = 0)

#Mutating the ORFtype columsns now together into a single column.
confdf %<>% mutate(nseqs = paste(cmp, `3p`, `5p`, int, sep = ", "))
confdf %<>% select(-c(cmp, `3p`, `5p`, int))


#Making protcat a factor column because I want the proteins grouped by
#the arm of the circadian cycle they belong to.
confdf$protcat <- factor(confdf$protcat, levels = c("CLK", "CYC", "TIM", "PER", 
                                                      "CRY1", "CRY2", "PDP1e", "VRI", 
                                                      "REV-ERBa", "RORa"))
confdf %<>% arrange(protcat)

#Fixing the sample names.
confdf %<>% 
  mutate(spmatch = str_replace(spmatch, "grp$", "")) %>%
  mutate(spmatch = str_replace(spmatch, "_", " ")) %>%
  mutate(spmatch = str_replace(spmatch, "sp$", "sp."))


#Pivoting into final table form.
#Form with the counts separated like this.
conftab1 <- confdf %>% 
  pivot_wider(names_from = protcat, values_from = nseqs, values_fill = NA)
#Form with counts as just numbers.
conftab2 <- confdf %>%
  rowwise() %>%
  mutate(nseqs = sum(as.numeric(unlist(str_split(nseqs, ", "))))) %>%
  data.frame() %>%
  pivot_wider(names_from = protcat, values_from = nseqs, values_fill = NA)
#Form as presence/absence matrix.
conftab3 <- confdf %>%
  rowwise() %>%
  mutate(nseqs = ifelse(sum(as.numeric(unlist(str_split(nseqs, ", ")))) > 0, "Yes", NA)) %>%
  data.frame() %>%
  pivot_wider(names_from = protcat, values_from = nseqs, values_fill = NA)



#Creating table.
#spmatch column name will be renamed on the fly
conftab_final <- conftab2

#Adding in phylum informationa and order.
conftab_final %<>% 
  mutate(disp_org = spmatch) %>% 
  full_join(taxdat, by = "disp_org") %>%
  arrange(ord) %>%
  select(-c(ord, disp_org))

conftab_final %<>% select(spmatch, disp_phy, !starts_with("spmatch|disp_phy"))

#Footnote:
genft <- "Numbers indicate number of candidates found. NAs indicate no candidates found. Circadian clock protein abbreviations: CLK (Circadian locomotor output cycles kaput), CYC (Cycle), PER (period), TIM1 (Timeless) CRY1 (Cryptochrome 1), CRY2 (Cryptochrome 2), PDP1e (PAR Domain Protein 1 epsilon), VRI (Vrille), REV-ERBa (Nuclear Receptor Subfamily 1 Group D Member 1/NR1D1), RORa (RAR Related Orphan Receptor A; NR1F1/Nuclear receptor subfamily 1 group F member 1). Phyla: Ar - Arthropoda; An - Annelida; Ph - Phoronida; Ch - Chordata; Ec - Echinodermata; Cn - Cnidaria."
#Writing to object
conftab_out <- conftab_final %>%
  dplyr::rename("Organism" = spmatch, "Phylum" = disp_phy) %>%
  #arrange(Organism) %>%
  kbl(., format = "latex", booktabs = T, align = "lcccccccccccc", linesep = "") %>%
  column_spec(1, italic = TRUE) %>%
  add_header_above(c(" " = 2, "Arthropod" = 8, "Mammalian" = 2)) %>%
  row_spec(row = 0, italic = FALSE, bold = TRUE) %>%
  footnote(general = genft, threeparttable = T) #%>%
  landscape()
#Here's the core latex table.
conftab_out

#Using save_kable with keep_tex = TRUE as a proxy to write the table out.
#Does not actually create the PDF since the GWDG R installation is missing
#dependencies.
#But the tex file does get written. This is all that's needed.
cat(conftab_out, file = paste0(visoutdir, "/", "of_cands_table.tex"))
save_kable(conftab_out, 
           file = paste0(visoutdir, "/", "of_cands_table.pdf"), 
           keep_tex = TRUE, self_contained = TRUE)

rm(confdf, conftab1, conftab2, conftab3, conftab_final, genft, conftab_out)




#--------------------------------------------------------------------------------


#DOMAIN STRUCTURE VISUALIZATIONS


#For the domain structure visualizations I'll be using seqvisr::pdomvisr().

#The annotations for this will be pulled from ipscands and ipsrefs.

#All I really need to do is get the protcat groupings of references and candidates
#from ccdf (and also prepare the names and whatnot properly), and it's good to go.
visdf <- ccdf %>% 
  select(c(protcat, refseq, matchseq, spref, spmatch, refmatcat))

#I need to pivot this into something longer so that references and matches are categories
#in columns, and not their own columns as they are right now.
#I need a custom spec for this.
myspec <- tribble(
  ~.name, ~.value, ~seq_class,
  "refseq", "seq", "ref",
  "matchseq", "seq", "match",
  "spref", "sp", "ref",
  "spmatch", "sp", "match"
)
#Pivotting
visdf %<>% pivot_longer_spec(spec = myspec)
rm(myspec)

#Just checking one last time to make sure none of the
#matches are repeated.
cat("Checking if there are any matches that occur more than once!!\n")
visdf %>% 
  filter(!str_detect(seq, "^SOIREF_|^OTHREF_")) %>% 
  group_by(sp) %>% 
  count(seq) %>% 
  ungroup() %>% 
  filter(n > 1) %>% {if(nrow(.) > 0) { cat("Yes")} else { cat("No\n") }}



#Also need to filter out the non-SOIREF__TRINITY rows at this point.
visdf %<>% filter(refmatcat == "SOIREF__TRINITY")

#Now applying distinct() to visdf to clear out the duplicated references.
visdf %<>% 
  group_by(sp) %>%
  distinct(seq, .keep_all = TRUE) %>%
  ungroup()

#I need to retain the original sequence identifiers + the species names
#as a unique key for merging in the InterProScan data, so doing that first in
#both cases.
visdf %<>% mutate(mergecol = paste0(sp, "__", seq))



#Munging the seq column to clean up the sequence identifiers.
#Its easier to deal with the references and matches as separate
#data.frames.
visrefs <- visdf %>% filter(str_detect(seq, "^SOIREF_|^OTHREF_"))
vismats <- visdf %>% filter(!str_detect(seq, "^SOIREF_|^OTHREF_"))


#Munging the references first.
visrefs %<>% mutate(seq = str_extract(seq, "[A-Za-z0-9]+(?=_[a-z]+$)"))
#visrefs %<>% mutate(seq = paste0("REF ", str_replace(sp, "_", " "), " ", seq))
visrefs %<>% mutate(seq = paste0(str_replace(sp, "_", " "), " ", seq, " REF"))

#Munging matches.
#Removing the length identifiers from the sequence name.
vismats %<>% mutate(seq = str_replace(seq, "_len_.*$", ""))
#Mutating the ORF information into single letter identifiers.
vismats %<>%
  mutate(seq = case_when(
    str_detect(seq, "_complete") ~ str_replace(seq, "_ORFtype.*", " (cmp)"),
    str_detect(seq, "_internal") ~ str_replace(seq, "_ORFtype.*", " (int)"),
    str_detect(seq, "_3prime") ~ str_replace(seq, "_ORFtype.*", " (3p)"),
    str_detect(seq, "_5prime") ~ str_replace(seq, "_ORFtype.*", " (5p)")
  ))
#Adding in the sample names
vismats %<>% 
  mutate(seq = paste0(str_replace(sp, "_", " "), " ", seq))

#Binding the two data.frames back together.
visdf <- bind_rows(visrefs, vismats)
rm(visrefs, vismats)


#--------------------------------------------------------------------------------


#Preparing ipsrefs and ipscands
#I can actually row bind these into a single data.frame.
ipsdat <- bind_rows(ipsrefs, ipscands)

#Munging filename and prot_acc cols to create mergecol.
ipsdat %<>% 
  mutate(filename = str_extract(basename(filename), "^[A-Za-z]+_[a-z]+")) %>%
  mutate(mergecol = paste0(filename, "__", prot_acc))

#Retaining columns.
ipsdat %<>% select(c(mergecol, seq_len, analysis_plat, signature_acc, 
                     signature_desc, start_loc, stop_loc))

rm(ipsrefs, ipscands)



#Also need to add the sites data to this.
ipsdat_sites <- bind_rows(ipsrefs_sites, ipscands_sites)

#Munging filename and prot_acc cols to create mergecol.
ipsdat_sites %<>% 
  mutate(filename = str_extract(basename(filename), "^[A-Za-z]+_[a-z]+")) %>%
  mutate(mergecol = paste0(filename, "__", prot_acc))

#Retaining columns.
ipsdat_sites %<>% select(c(mergecol, seq_len, analysis_plat, signature_acc, 
                     signature_desc, start_loc, stop_loc))

rm(ipsrefs_sites, ipscands_sites)

#Row-binding ipsdat_sites with ipsdat.
#I want to add a column to distinguish the sites from the domains.
ipsdat$annot_type <- "dom"
ipsdat_sites$annot_type <- "site"
ipsdat <- bind_rows(ipsdat, ipsdat_sites)
rm(ipsdat_sites)

#Replace the "grp" in the names in ipsdat before passing on to merge.
#This is because it's been replaced in visdf already.
ipsdat %<>% mutate(mergecol = str_replace(mergecol, "grp__", "__"))


#--------------------------------------------------------------------------------

#flps, nls if needed.

#--------------------------------------------------------------------------------
#Merging.
visdf <- merge(visdf, ipsdat, by = "mergecol", all.x = TRUE)
rm(ipsdat)

#There are some repeated annotations in here still.
#E.g., Drosophila CYC's CDD annotations are repeated.
#Filtering these out.
visdf %<>% 
  group_by(mergecol, annot_type, analysis_plat, signature_acc) %>% 
  distinct(start_loc, .keep_all = TRUE) %>%
  ungroup()

visdf <- data.frame(visdf, stringsAsFactors = FALSE)


#--------------------------------------------------------------------------------

#PLOTTING PLANNING SECTION.

#FIRST IMPORTANT THING.
#Are there any sequences with no annotations?
visdf %>% filter(is.na(signature_acc))
#No, but if there were, the annotation columns would
#have to be set up appropriately in order for seqvisr()
#to be able to handle them.


#First I will take a look at how the protein categories are
#covered by the various annotation platforms.
visdf %>% 
  group_by(analysis_plat) %>% 
  count(protcat) %>% 
  pivot_wider(names_from = protcat, values_from = n)


#Only need Pfam + CDD-3.18 for annotations.
#I'll need CDD or SMART to fill in for CLK alone,
#as it's the only case where Pfam messes up a reference.
#I'll keep CDD.
visdf %<>% 
  filter(!(analysis_plat %in% c("ProSitePatterns", "Phobius", "SMART",
                                "Gene3D", "ProSiteProfiles", "PANTHER")))

#--------------------------------------------------------------------------------


#Since I have both Pfam and CDD annotations, I need to 
#"merge" the annotations in such a way that the CDD annotations
#get "painted in" only when an equivalent Pfam annotation is missing.

#For this, I first need to split away the sites annotations, as they
#would just get subsumed otherwise.
visdf_sites <- visdf %>% filter(annot_type == "site")
visdf %<>% filter(annot_type != "site")


#Merging overlapping domain sets.
#From https://stackoverflow.com/q/15235821
for(i in 1:length(unique(visdf$mergecol))){
  
  #i <- 1
  curdat <- visdf %>% filter(mergecol == unique(visdf$mergecol)[i])
  
  #To merge domain sets.
  #From https://stackoverflow.com/q/15235821
  ir <- IRanges(curdat$start_loc, curdat$stop_loc)
  #Assigning overlapping domains to the same domgrp set.
  curdat$domgrp <- S4Vectors::subjectHits(IRanges::findOverlaps(ir, IRanges::reduce(ir)))
  
  #Calculating the lengths of the domains
  curdat$domlen <- curdat$stop_loc - curdat$start_loc
  
  #For each group of domains, I want a consistent description.
  #So making sure it's the Pfam description (if available).
  #Doing this with the analysis_ord mutate + subsequent arrange() below
  #The "new" annotation will go into the sig_new column.
  curdat %<>%
    group_by(domgrp) %>%
    mutate(analysis_ord = ifelse(analysis_plat == "Pfam", 1, 2)) %>%
    arrange(analysis_ord, .by_group = TRUE) %>%
    mutate(sig_new = signature_desc[1]) %>%
    #slice_max(domlen, n = 1, with_ties = FALSE) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  if(i == 1){
    outdf <- curdat
  } else{
    outdf <- bind_rows(outdf, curdat)
  }
  
}
rm(i, ir, curdat)

#How many protein annotation sources do we have?
outdf %>% 
  group_by(protcat) %>% 
  distinct(analysis_plat) %>% 
  arrange(analysis_plat, .by_group =  TRUE) %>%
  mutate(plats = paste0(analysis_plat, collapse = ",")) %>%
  ungroup() %>% 
  group_by(protcat) %>%
  distinct(plats)
#Looks good.

#Overwriting visdf with outdf now.
bkup <- visdf
visdf <- outdf
#rm(outdf)

#In case I ran slice_max(domlen), I might have to run this
#code below.
#So I will filter out CDD and ProSiteProfiles also now.
#visdf %<>% filter(!(analysis_plat %in% c("CDD", "ProSiteProfiles")))



#--------------------------------------------------------------------------------

#Need to clean up the signature descriptions before plotting.
visdf %>% distinct(sig_new)

#Only thing that really needs to be cleaned up is the multiple PAS
#domain annotations.
visdf %<>% 
  mutate(sig_new = ifelse(str_detect(sig_new, "^PAS"), "PAS domain", sig_new))

#Removing commas from the description.
visdf %<>% mutate(sig_new = str_replace_all(sig_new, ",", ""))

#Checking.
visdf %>% distinct(sig_new)

#Need to come up with succinct names for these (which I can then elaborate in the legend)
doms <- visdf %>% distinct(sig_new)

doms %<>% mutate(sig_small = case_when(
  sig_new == "Zinc finger C4 type (two domains)" ~ "ZF",
  sig_new == "Ligand-binding domain of nuclear hormone receptor" ~ "LBD",
  sig_new == "Helix-loop-helix DNA-binding domain" ~ "bHLH",
  sig_new == "PAS domain" ~ "PAS",
  sig_new == "DNA photolyase" ~ "DNA-p",
  sig_new == "FAD binding domain of DNA photolyase" ~ "FAD-b",
  sig_new == "Basic region leucine zipper" ~ "bZIP",
  sig_new == "Period protein 2/3C-terminal region" ~ "Per-C",
  sig_new == "Timeless protein" ~ "Tim-N",
  sig_new == "Timeless PAB domain" ~ "PAB",
))

visdf <- merge(visdf, doms, by = "sig_new", all.x = TRUE)

#--------------------------------------------------------------------------------

#Assigning colors for the domains.
#I don't need 10 unique colors despite having 10 unique domains.
#Just going to recycle colors, but in such a manner that the colors
#assigned to co-occurring domains are distinct.

#Finding out which domains co-occur together.
# domcols <- visdf %>% 
#   group_by(protcat) %>%
#   distinct(sig_small) %>%
#   mutate(domgrps = paste0(sig_small, collapse = ",")) %>%
#   ungroup() %>%
#   distinct(sig_small, .keep_all = TRUE) %>%
#   select(-protcat)
# domcols

domcols <- visdf %>% 
  group_by(protcat) %>%
  distinct(sig_small) %>%
  mutate(domgrps = paste0(sig_small, collapse = ",")) %>%
  ungroup() %>%
  separate_rows(domgrps, sep = ",") %>%
  group_by(sig_small) %>%
  mutate(domgrps = paste0(unique(domgrps), collapse = ",")) %>%
  distinct(sig_small, .keep_all = TRUE) %>%
  ungroup() %>%
  select(-protcat)

domcols

#Orange: #f58231
#I will just assign colors manually.
domcols %<>% 
  mutate(featcols = case_when(
    sig_small == "bZIP" ~ "maroon",
    sig_small == "DNA-p" ~ "#e6194B",
    sig_small == "FAD-b" ~ "#4363d8",
    sig_small == "bHLH" ~ "#4363d8",
    sig_small == "LBD" ~ "#40B0A6",
    sig_small == "PAS" ~ "seagreen",
    sig_small == "Per-C" ~ "#f58231", 
    sig_small == "PAB" ~ "#3cb44b",
    sig_small == "Tim-N" ~ "#4363d8",
    sig_small == "ZF" ~ "orchid4"
  ))

#Selecting columns for merging in the colors.
domcols %<>% select(c(sig_small, featcols))

#Merging in the colors.
#WILL USE A NEW DATA.FRAME CALLED VISOUT.
visout <- merge(visdf, domcols, by = "sig_small", all.x = TRUE)
visout <- data.frame(visout, stringsAsFactors = FALSE)

#These colors are acceptable (from having tested their combinations manually).

#Adding in heights and offsets.
#Only really necessary if I'm also plotting the sites,
#but will keep these in anyway.
visout %<>% mutate(feath = ifelse(annot_type == "site", 0.4, 0.6),
       offset = 0)



#--------------------------------------------------------------------------------


#PLOTTING SECTION.


#I want to add the InterProScan feature accession to the
#feature signature descriptor columns (both sig_new and sig_small).
visout %<>% mutate(sig_new = paste0(sig_small, " (", sig_new, " ", signature_acc, ")"))

#Sigset is a data.frame tying all unique sig_small's to their sig_new's.
#This is so that I can use this data.frame to create legend labels that
#are sig_small's but the in-place feature descriptors are sig_new's.
sigset <- visout %>% 
  select(sig_small, sig_new, featcols) %>%
  distinct(sig_small, .keep_all = TRUE)



#Selecting plotting relevant columns.
visout %<>% 
  select(c(protcat, analysis_plat, seq, seq_len, offset, 
           feath, featcols, sig_small, start_loc, stop_loc))




#----------------------------------------------------------------


##IMAGE PLOTTING FINAL##



#Using a for loop to plot the images.
for(i in 1:length(unique(visout$protcat))){
  
  #      i <- 1
  #selprot <- "PER"
  selprot <- unique(visout$protcat)[i]
  
  cat("Plotting ", selprot, "!!\n")
  
  #Selecting the plot data.
  pltdat <- visout %>%
    filter(protcat == selprot & (analysis_plat %in% c("Pfam", "CDD"))) %>%
    #arrange(ord) %>%
    select(-c(protcat, analysis_plat))
  
  #In order to have sig_small as the in-feature label, and sig_new as the
  #corresponding legend labels, I am going to have to intercept the output
  #of seqvisr::tsvtogginp_multi(), and set these values in the appropriate
  #columns. I can then pass this to the pdomvisr() plotting function.
  
  #So creating a data.frame for the intermediate data.
  pltdat <- seqvisr::tsvtogginp_multi(pltdat)
  
  #The label column already contains the right value (i.e., the sig_small value)
  #I just need to update the posdesc values.
  pltdat %<>% mutate(sig_small = posdesc)
  pltdat %<>% left_join(y = sigset, by = "sig_small")
  pltdat %<>% mutate(posdesc = sig_new) %>% select(-c(sig_small, sig_new, featcols))
  
  
  #Setting up italics in the sequence names, these will be interpreted by
  #ggtext::element_markdown() in theme().
  pltdat %<>% mutate(prot_acc = str_replace(prot_acc, "(^[A-Za-z]+\\s)([a-z]+)(.*$)", "\\1\\2*\\3"),
                     prot_acc = paste0("*", prot_acc))
  
  pltdat %<>% 
    mutate(seq_class = ifelse(str_detect(prot_acc, "REF"), 0, 1)) %>% 
    arrange(seq_class) %>% 
    mutate(prot_acc = factor(prot_acc, levels = rev(unique(prot_acc)))) %>%
    select(-seq_class)
  
  #I want to str_wrap posdesc as some of these strings can be really long.
  pltdat %<>% mutate(posdesc = str_wrap(posdesc, width = 30))
  
  #Adding the . after sp in prot_acc.
  pltdat %<>% mutate(prot_acc = str_replace(prot_acc, "(?<=\\s)(sp)(?=\\*)", "\\1\\."))
  
  #Adding in the phylum data and ordering in the taxonomic order.
  pltdat %<>% 
    mutate(disp_org = str_extract(prot_acc, "(?<=\\*)[A-Za-z]+\\s[a-z\\.]+(?=\\*)")) %>%
    full_join(taxdat, by = "disp_org") %>%
    mutate(ord = ifelse(is.na(ord), 0, ord)) %>%
    mutate(prot_acc = ifelse(!str_detect(prot_acc, "\\sREF"), paste0(prot_acc, " (", disp_phy, ")"), prot_acc))
  pltdat %<>% arrange(ord)
  #This stage can add NA rows since not all taxa are represented in all plots.
  #Filtering out NA rows.
  pltdat %<>% filter(!is.na(prot_acc))
  pltdat %<>% select(-c(disp_phy, disp_org, ord))
  
  #Factoring the prot_acc column and reverse ordering it.
  #pltdat %<>% mutate(ordloc = row_number())
  pltdat$prot_acc <- forcats::fct_inorder(pltdat$prot_acc)
  pltdat$prot_acc <- forcats::fct_rev(pltdat$prot_acc)
  #pltdat$prot_acc <- factor(pltdat$prot_acc, levels = pltdat$ordloc)
  #pltdat %<>% select(-ordloc)
  
  #hbase = 0.2, label_size = 4
  plt <- pdomvisr(inpdat = pltdat, label_size = 5, hbase = 0.4, ylabel = "") + #size = 2, hbase = 0.2
    theme(text = ggplot2::element_text(size = 18, face = "bold"), #16
          legend.position = "top",
          legend.direction = "horizontal",
          legend.justification = "left",
          legend.key.size = unit(0.5, "cm"), 
          legend.text = ggplot2::element_text(size = 16, face = "bold", hjust = 0), #size = 14
          legend.title = ggplot2::element_text(size = 24, face = "bold"), #size = 20
          axis.text.y = ggtext::element_markdown(size = 20), #size = 18
          legend.box.just = "left", 
          legend.text.align = 0,
          plot.title.position = "panel", 
          plot.title = element_text(vjust = 0.5)) + 
    labs(title = "") +
    guides(fill = guide_legend(title = paste0(selprot, " features"), nrow = 1, hjust = 0))
  plt
  
  
  plt$plot_env$inpdf %<>% 
    mutate(disp_org = str_extract(prot_acc, "(?<=\\*)[A-Za-z]+\\s[a-z\\.]+(?=\\*)")) %>%
    full_join(taxdat, by = "disp_org") %>%
    mutate(ord = ifelse(is.na(ord), 0, ord)) %>%
    mutate(prot_acc = ifelse(!str_detect(prot_acc, "\\sREF"), paste0(prot_acc, " (", disp_phy, ")"), prot_acc)) %>%
    arrange(ord) %>%
    select(-c(disp_phy, disp_org, ord))
  
  plt
  
  #Saving PNG.
  outfile <- paste0("seqvisr_", selprot, "_all.png")
  outfile <- paste0(visoutdir, "/", outfile)
  ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "png")
  #Saving PDF.
  outfile <- paste0("seqvisr_", selprot, "_all.pdf")
  outfile <- paste0(visoutdir, "/", outfile)
  ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = cairo_pdf)
  #Saving SVG.
  outfile <- paste0("seqvisr_", selprot, "_all.svg")
  outfile <- paste0(visoutdir, "/", outfile)
  ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "svg")
  
  cat("Done plotting ", selprot, "!!\n")
  
}






#--------------------------------------------------------------------------------
