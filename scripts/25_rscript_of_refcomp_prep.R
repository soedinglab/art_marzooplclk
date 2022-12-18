#This script does some preparatory work to compare
#the existing OrthoFinder output with three reference species
#to the output from a run with more reference species to see
#if this affects detection rates.

#Step 1:
#This script uses a pre-generated API URL from UniProt
#to download a TSV file of all non-redundant metazoan
#reference proteomes from UniProt except those of D. melanogaster,
#D. plexippus, and M. musculus.
#The FASTA files for these proteomes have already been downloaded
#and are in local storage.

rm(list = ls())

library(ape)
library(taxize)
library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)

mypath <- "/path/to/of_refcomp"
setwd(mypath)


#Downloading additional reference proteomes identifiers from UniProt
updl_url <- "https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd%2Clineage&format=tsv&query=%28%28taxonomy_id%3A33208%29%20AND%20%28busco%3A%5B95%20TO%20%2A%5D%29%20NOT%20%28taxonomy_id%3A7227%29%20NOT%20%28taxonomy_id%3A10090%29%20NOT%20%28taxonomy_id%3A13037%29%20NOT%20%28proteome_type%3A2%29%20NOT%20%28proteome_type%3A3%29%20NOT%20%28proteome_type%3A4%29%29"
dlfile <- paste0(mypath, "/", "up_metazoans_list.tsv")
download.file(url = updl_url, destfile = dlfile, method = "curl", cacheOK = TRUE)
updat <- data.table::fread(updl_url, sep = "\t", header = TRUE)

#Setting this up for merging.
updat %<>% 
  mutate(samplename = str_extract(Organism, "^[A-Za-z]+\\s[a-z]+")) %>%
  mutate(samplename = str_replace(samplename, "\\s+", "_")) %>%
  mutate(samplename = paste0(samplename, "_", `Proteome Id`, "_", `Organism Id`)) %>%
  rename("ncbi_taxid" = `Organism Id`) %>%
  select(c(samplename, ncbi_taxid))

#Appending a "_nonrefs" string to these filenames to make them easily identifiable.
#This has already been added to the actual filenames also.
updat %<>% mutate(samplename = paste0(samplename, "_nonrefs"))

#Reading in table of species names and NCBI taxonomy IDs I prepared manually.
taxdat <- read.table(paste0(mypath, "/", "orthofinder_ncbi_taxonomy_for_correct_species_tree.csv"), sep = ",", header = TRUE)
taxdat %<>% select(c(samplename, ncbi_taxid))

#Munging
#taxdat %<>% select(c(organism, ncbi_taxid))
taxdat %<>% mutate(ncbi_taxid = as.numeric(str_extract(ncbi_taxid, "\\d+")))


#Merging all identifiers.
taxdat <- bind_rows(taxdat, updat)
rm(updat)

#Fetching taxonomy tree from NCBI
ncbidat <- classification(taxdat$ncbi_taxid, db = "ncbi")
ncbidat <- class2tree(ncbidat, check = TRUE)

plot(ncbidat)
ncbidat$phylo$tip.label


for(i in 1:length(ncbidat$names)){
  if(taxdat$ncbi_taxid[i] == ncbidat$names[i]){
    ncbidat$phylo$tip.label[i] <- taxdat$samplename[i]
  }
}

plot(ncbidat)

write.tree(ncbidat$phylo, file = paste0(mypath, "/", "of_sp_tree_extrefs.nwk"))



