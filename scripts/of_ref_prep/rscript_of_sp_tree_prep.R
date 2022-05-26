#OrthoFinder needs a species phylogenetic tree to increase its output accuracy.
#NCBI can only give a taxonomy tree.
#So following this: https://taylorreiter.github.io/2017-07-28-Taxonomy-from-Species-Name-in-R/
#

library(ape)
library(taxize)
library(magrittr)
library(tidyr)
library(stringr)
library(dplyr)

mypath <- "/path/to/of_ref_prep" #SET THIS PATH!!

#Reading in table of species names and NCBI taxonomy IDs I prepared manually.
taxdat <- read.table(paste0(mypath, "/", "orthofinder_ncbi_taxonomy_for_correct_species_tree.csv"), sep = ",", header = TRUE)

#Munging
#taxdat %<>% select(c(organism, ncbi_taxid))
taxdat %<>% mutate(ncbi_taxid = as.numeric(str_extract(ncbi_taxid, "\\d+")))

#Will not retain Synechococcus and Neurospora.
#taxdat %<>% filter(!str_detect(organism, "Synechococcus|Neurospora"))

ncbidat <- classification(taxdat$ncbi_taxid, db = "ncbi")
ncbidat <- class2tree(ncbidat, check = TRUE)

plot(ncbidat)

ncbidat$phylo$tip.label


#tst <- ncbidat


for(i in 1:length(ncbidat$names)){
  if(taxdat$ncbi_taxid[i] == ncbidat$names[i]){
    ncbidat$phylo$tip.label[i] <- taxdat$samplename[i]
  }
}

plot(ncbidat)

write.tree(ncbidat$phylo, file = paste0(mypath, "/", "of_sp_tree.nwk"))
