#Script to check how many of our organisms have
#extant transcriptomes (or genomes) in the NCBI databases.

rm(list = ls())

library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)
library(rentrez)


#Reading in file containing species names and NCBI taxids.
mypath <- "/home/owner/Nextcloud/laptop_rplace/eichele/newrun/scripts/ncbi_assem_counts/species_taxids.csv"
spdat <- read.table(file = mypath, sep = ",", header = TRUE)

#Appending "txid" to the taxids since this is missing.
spdat %<>% mutate(ncbi_taxid = paste0("txid", ncbi_taxid))
#spdat %<>% mutate(ncbi_taxid = ifelse(!str_detect(organism, "sp\\.$"), ncbi_taxid, NA))

#Getting the transcriptome and genome assembly
#counts from NCBI.
for(i in 1:nrow(spdat)){
  
  sterm <- paste0("tsa-master[PROP] AND ", spdat$ncbi_taxid[i], "[ORGN]")
  tmp <- entrez_search(db = "nuccore", term = sterm)
  spdat$num_transcriptomes[i] <- tmp$count
  
  cat(spdat$organism[i], " Transcriptomes: ", tmp$count, "\n")
  
  sterm <- paste0(spdat$ncbi_taxid[i], "[ORGN]")
  tmp <- entrez_search(db = "genome", term = sterm)
  spdat$num_genomes[i] <- tmp$count
  
  cat(spdat$organism[i], " Genomes: ", tmp$count, "\n")
  
}

#Number of assemblies, genomic OR transcriptomic
spdat %<>% mutate(num_tot_assem = num_transcriptomes + num_genomes)

#How many have any assemblies at all?
spdat %>% filter(num_tot_assem > 0) %>% nrow()

#Writing this out.
write.table(spdat, file = "ncbi_assem_counts.csv", sep = ",", row.names = FALSE)
