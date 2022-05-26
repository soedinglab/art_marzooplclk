#Script to prepare a species tree figure for the organisms included
#in this analysis.


rm(list = ls())

#Fonts.
library("extrafont")
#See https://stackoverflow.com/a/68642855/9494044
#library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#Run font_import for the first time after package installation only!!
#extrafont::font_import()

#For exporting EPS images mainly.
loadfonts(device = "postscript")

#For getting high quality PNGs from PDFs (ggsave problems).
library(pdftools)

#Tree-related packages.
library(ggplotify)
library(taxize)
library(treeio)
library(ape)
library(tidytree)
library(ggtree)

#General-purpose.
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)


#Main path.
mypath <- "/home/owner/Nextcloud/laptop_rplace/eichele/newrun/rscript_figstabs"

setwd(mypath)

#----------------------------------------------------------------

##SPECIES TREE FIGURE##

#Reading in table of species names and NCBI taxonomy IDs I prepared manually.
taxdat <- read.table(paste0(mypath, "/", "species_data.csv"), sep = ",", header = TRUE)

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

#Renaming the tips in the tree object.
for(i in 1:length(ncbitree$names)){
  for(j in 1:length(taxdat$ncbi_taxid)){
    if(taxdat$ncbi_taxid[j] == ncbitree$names[i]){
      ncbitree$phylo$tip.label[i] <- taxdat$organism[j]
    }
  }
}
rm(i, j)
#Extracting the tree-relevant phylo object from ncbidat.
ncbitree <- ncbitree$phylo


#----
#Plotting species tree.

#Now I will use ncbitree and taxdat via ggtree() to plot a phylogenetic tree
#with colored tip labels.
#First column NEEDs to be a column corresponding to the tip labels.
#So rearranging taxdat.
taxdat %<>% select(names(taxdat)[2], names(taxdat)[1], names(taxdat)[-c(1:2)])



#Need to change the ; in dataset_indiv to ,.
#taxdat %<>% mutate(dataset_indiv = str_replace(dataset_indiv, ";", ", "))

#Creating new column for the plotlabels.
#I want to add a * to the plabel if the sample was processed in 2020.
taxdat %<>% 
  mutate(plabel = ifelse(str_detect(sampling_time, "2020"), paste0(organism, "*"), organism))

#Adding in the dataset data to this label.
taxdat %<>%
  mutate(plabel = paste0(plabel, " (", dataset_indiv, ")"))

#Adding padded labels.
#taxdat %<>% mutate(plabel = label_pad(organism, pad = "."))
#taxdat$plabel <- str_pad(taxdat$organism, max(nchar(taxdat$organism))+10, side = "left", pad = "·")
#Also need to pad the datasets column.
#taxdat$dataset_indiv <- str_pad(taxdat$dataset_indiv, max(nchar(taxdat$dataset_indiv)), side = "right", pad = " ")
taxdat %<>% mutate(plabel = str_pad(plabel, max(nchar(plabel))+10, side = "left", pad = "."))


#Generating basic plot.
treeplt <- ggtree(ncbitree) %<+% taxdat

#Fleshed-out plot.
#Excluding nodes #Exclude nodes 31, 30, 26, 29
#I found these nodes by first setting label = node in geom_label2
#and not subsetting anything except isTip.
filtnodes <- c(31, 30, 26, 29)


treeplt0 <- treeplt + 
  geom_tiplab(size = 6, fontface = 3, align = TRUE, mapping = aes(label = plabel), family='mono') +
  geom_tippoint(aes(shape = phylum, color = phylum), size = 6) + 
  geom_label2(aes(subset = !isTip & !(node %in% filtnodes), label = label), fill = "orange", size = 6, alpha = 0.8) + 
  xlim(-4, 100) +
  theme(legend.position = c(0.120, 0.775), 
        legend.key.size = unit(1, 'cm'), 
        legend.title = ggplot2::element_text(size = 20, face = "bold"), 
        legend.text = ggplot2::element_text(size = 15, face = "bold.italic"),
        legend.background = ggplot2::element_rect(fill = "gray80",
                                         size = 0.5, linetype = "solid", 
                                         colour = "black")) + 
  ggplot2::guides(color = guide_legend(title = "Phylum"), 
                  shape = guide_legend(title = "Phylum"))

treeplt0

#Writing out plot to file.
#PDF
treefile <- paste0(mypath, "/", "species_tree.pdf")
ggsave(filename = treefile, plot = treeplt0, dpi = 1200, width = 40, height = 20, units = "cm", limitsize = FALSE)
#PNG
#pdftools::pdf_convert(pdf = treefile, format = "png", dpi = 1200)
#EPS
#treefile <- paste0(mypath, "/", "species_tree.eps")
#ggsave(filename = treefile, plot = treeplt0, dpi = 1200, width = 40, height = 20, units = "cm", device = cairo_ps)


#----------------------------------------------------------------

# #For this version, if I am to use it, I should reset the plabel to just
# #point to the organism name, and just pad that.
# #Extended plot.
# 
# treeplt1 <- treeplt + 
#   geom_tiplab(size = 6, fontface = 3, align = TRUE, mapping = aes(label = plabel)) +
#   geom_tippoint(aes(shape = phylum, color = phylum), size = 6) + 
#   geom_label2(aes(subset = !isTip & !(node %in% filtnodes), label = label), fill = "orange", size = 6, alpha = 0.8) + 
#   xlim(-4, 80) +
#   ggplot2::guides(color = guide_legend(title = "Phylum"), 
#                   shape = guide_legend(title = "Phylum"))
# 
# treeplt1
# 
# #Adding in data columns to display as geom_facets.
# #ggtree is okay with this only in a very specific format,
# #in that each column needs to be its own data.frame with a
# #column called ID corresponding to the tip labels.
# 
# #Sampling time.
# facetdat1 <- taxdat %>% 
#   rename(id = organism, value = sampling_time) %>% 
#   select(c(id, value))
# treeplt2 <- facet_plot(treeplt1, mapping = aes(x = 0, label = value),
#                        hjust = -0.2,
#                        data = facetdat1, 
#                        geom = geom_text, 
#                        panel = "Sampling date")
# treeplt2
# 
# #Number of individuals.
# facetdat2 <- taxdat %>%
#   rename(id = organism, value = dataset_indiv) %>% 
#   select(c(id, value))
# treeplt3 <- facet_plot(treeplt2, mapping = aes(x = 0, label = value),
#                        hjust = -0.305,
#                        data = facetdat2, 
#                        geom = geom_text, 
#                        panel = "Datasets (individuals)")
# 
# treeplt3
# 
# #Need to move the legend into position BEFORE facet_widths.
# treeplt4 <- treeplt3 + theme(legend.position = c(0.10, 0.76))
# 
# #treeplt4
# #Setting facet label sizes here also.
# treeplt4 <- treeplt4 + 
#   theme(#legend.position = c(0.26, 0.78), 
#   legend.key.size = unit(1, "cm"),
#   legend.title = ggplot2::element_text(size = 20, face = "bold"), 
#   legend.text = ggplot2::element_text(size = 14, face = "bold.italic"),
#   legend.background = ggplot2::element_rect(fill = "gray80",
#                                             size = 0.5, linetype = "solid", 
#                                             colour = "black"), 
#   strip.text = ggplot2::element_text(size = 12, face = "bold")) 
# 
# 
# #, text = ggplot2::element_text(family = "Helvetica")
# 
# 
# 
# #Updating the facet widths to make everything size-meaningful.
# treeplt5 <- facet_widths(p = treeplt4, 
#              widths = c(`Tree` = 1.0, 
#                         `Sampling date` = 0.15, 
#                         `Datasets (individuals)` = 0.20))
# #Final plot.
# treeplt5
# 
# 
# 
# #Writing out plot to file.
# #For this extended figure, for whatever reason, the tiff and png versions do not render the
# #legend box properly (some words overflow the box).
# #So going to write out a PDF, and then use the pdftools package to render a high quality
# #PNG from that.
# #treefile <- paste0(mypath, "/", "species_tree_samplingdata.png")
# #ggsave(filename = treefile, plot = treeplt5, dpi = 600, width = 40, height = 20, units = "cm", limitsize = FALSE)
# 
# #PDF
# treefile <- paste0(mypath, "/", "species_tree_samplingdata.pdf")
# ggsave(filename = treefile, plot = treeplt5, dpi = 600, width = 40, height = 20, units = "cm", limitsize = FALSE)
# #PNG
# pdftools::pdf_convert(pdf = treefile, format = "png", dpi = 600)
# #EPS
# treefile <- paste0(mypath, "/", "species_tree_samplingdata.eps")
# ggsave(filename = treefile, plot = treeplt5, dpi = 600, width = 40, height = 20, units = "cm", device = cairo_ps)

#----------------------------------------------------------------

