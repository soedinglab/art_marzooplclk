#Script to check how many of our organisms have
#extant transcriptomes (or genomes) in the NCBI databases.

rm(list = ls())

library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)
library(rentrez)
#library(reutils)


#Reading in file containing species names and NCBI taxids.
mypath <- "/home/owner/Nextcloud/laptop_rplace/eichele/art_marzooplclk/scripts/ncbi_assem_counts/species_taxids.csv"
spdat <- read.table(file = mypath, sep = ",", header = TRUE)

#Appending "txid" to the taxids since this is missing.
spdat %<>% mutate(ncbi_txid = paste0("txid", ncbi_txid))
#spdat %<>% mutate(ncbi_txid = ifelse(!str_detect(organism, "sp\\.$"), ncbi_txid, NA))

#For the species that haven't been identified, inserting the genus txid.
# extdat <- read.table(text = "organism\tncbi_txid\nCorystes sp.\ttxid557244\nHyperia sp.\ttxid371513\nPoecilochaetus sp.\ttxid273030\n")
# extdat <- data.frame(organism = c("Corystes sp.", "Hyperia sp.", "Poecilochaetus sp."), 
#                      ncbi_txid = c("txid557244", "txid371513", "txid273030"), 
#                      comment = c("used_genus_txid", "used_genus_txid", "used_genus_txid"))

#Merging and retaining only rows with NOT NA taxids.
# spdat %<>% bind_rows(extdat) %>% filter(!str_detect(ncbi_txid, "NA$"))


#Creating additional columns that'd be needed.
spdat$num_transcriptomes <- NA
spdat$num_genomes <- NA
spdat$num_refseq_protein <- NA
spdat$num_refseq_nucleotide <- NA
spdat$bioproject_accession_transcriptomes <- NA
spdat$bioproject_accession_genomes <- NA



#Getting the transcriptome and genome assembly
#counts from NCBI.
for(i in 1:nrow(spdat)){
  
  #  i <- 17
  #  i <- 4
  #  i <- 3
  
  #Getting counts of transcriptomes
  sterm <- paste0("tsa-master[PROP] AND ", spdat$ncbi_txid[i], "[ORGN]")
  tmp <- entrez_search(db = "nuccore", term = sterm)
  #spdat$num_transcriptomes[i] <- tmp$count
  #Getting associated bioproject accessions
  if(length(tmp$ids) > 0){
    tmp1 <- entrez_link(dbfrom = "nuccore", id = tmp$ids, db = "bioproject")
    if(all(class(tmp1) == c("elink", "list"))){
      biopids <- tmp1$links$nuccore_bioproject_tsamaster
    }
    if(all(class(tmp1) == c("elink_list", "list"))){
      biopids <- tmp1[[1]]$links$nuccore_bioproject
    }
    if(length(biopids) > 0){
      tmp2 <- entrez_fetch(db = "bioproject", id = biopids, rettype = "xml", parsed = TRUE)
      tmp3 <- entrez_summary(db = "bioproject", id = biopids)
      tmp4 <- extract_from_esummary(tmp3, "project_acc")
      tmp4 <- paste0(tmp4, collapse = ",")
      spdat$bioproject_accession_transcriptomes[i] <- tmp4
    } else {
      spdat$bioproject_accession_transcriptomes[i] <- NA
    }
  } else {
    spdat$bioproject_accession_transcriptomes[i] <- NA
  }
  
  #Counting this way because multiple nuccore records can point
  #to the same BioProject.
  spdat$num_transcriptomes[i] <- str_count(spdat$bioproject_accession_transcriptomes[i], "PRJ[A-Z0-9]+")
  if(is.na(spdat$num_transcriptomes[i])){
    spdat$num_transcriptomes[i] <- 0
  }
  
  cat(spdat$organism[i], " Transcriptomes: ", tmp$count, "\n")
  
  
  
  #Getting counts of genomes
  sterm <- paste0(spdat$ncbi_txid[i], "[ORGN]")
  tmp <- entrez_search(db = "genome", term = sterm)
  spdat$num_genomes[i] <- tmp$count
  #Getting associated bioproject accessions
  if(length(tmp$ids) > 0){
    tmp1 <- entrez_summary(db = "genome", id = tmp$ids)
    tmp1 <- paste0(tmp1$project_accession, collapse = ",")
    spdat$bioproject_accession_genomes[i] <- tmp1
  } else {
    spdat$bioproject_accession_genomes[i] <- NA
  }
  
  cat(spdat$organism[i], " Genomes: ", tmp$count, "\n")
  
  
  #Getting RefSeq protein counts.
  sterm <- paste0(spdat$ncbi_txid[i], "[ORGN] AND srcdb_refseq[property]")
  tmp <- entrez_search(db = "protein", term = sterm)
  spdat$num_refseq_protein[i] <- tmp$count
  
  cat(spdat$organism[i], " RefSeq proteins: ", tmp$count, "\n")
  
  
  #Getting RefSeq nucleotide sequence counts.
  sterm <- paste0(spdat$ncbi_txid[i], "[ORGN] AND srcdb_refseq[property]")
  tmp <- entrez_search(db = "nucleotide", term = sterm)
  spdat$num_refseq_nucleotide[i] <- tmp$count
  
  cat(spdat$organism[i], " RefSeq nucleotide sequences: ", tmp$count, "\n")

  
}

#Number of assemblies, genomic OR transcriptomic
spdat %<>% mutate(num_tot_assem = num_transcriptomes + num_genomes)

#How many have any assemblies at all?
spdat %>% filter(num_tot_assem > 0) %>% nrow()

#Writing this out.
write.table(spdat, file = "ncbi_assem_counts.csv", sep = ",", row.names = FALSE)


#-----------------------------------------------------------------------------------------------------#

