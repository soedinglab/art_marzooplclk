#Script to take all cases where multiple circadian clock candidates were found,
#align them to each other, and compare the results.
#Alignments will be done with the msa package, and visualizations with seqvisr.

##IMPORTANT##
#Make sure the FASTA files from the GitHub repository 
#under "circadian_clock_candidates/fas_by_type" are placed
#under mypath in "rscript_ccsel/fas_by_type"
#Or if rscript_ccsel's script has been run, they would automatically
#be here in this location.
####

rm(list = ls())

#Packages
#devtools::install_github("vragh/seqvisr")
library(seqvisr)
#BiocManager::install("msa")
library(msa)
library(Biostrings)
library(gtools) #For permutations and combinations
library(ggplot2)
library(magrittr)
library(stringr)
library(dplyr)
library(tidyr)


#Main paths
mypath <- "/path/to/outputs"


#Defining and creating the output directory (if it doesn't exist already)
outpath <- paste0(mypath, "/", "rscript_ccseqcomp")
if(!dir.exists(outpath)){
  dir.create(outpath)
}
setwd(outpath)

#Input files are the FASTA files written out by the 21b_rscript_ccsel_mainjob.R script.
#These files are under rscript_ccsel/fas_by_type.
cands <- paste0(mypath, "/", "rscript_ccsel/fas_by_type")
#Reading in the sequences
cands <- seqvisr::fasdirdf(path = cands, seqtype = "AA")

#The reference sequences aren't needed for this analysis.
#Dropping here already so that the species name and protein category extraction are
#made easier because of the consistent naming pattern available thereafter.
cands %<>% filter(!str_detect(seqname, "^SOIREF"))

#Extracting the clock protein category and species name.
cands %<>% 
  rowwise() %>% 
  mutate(species = unlist(str_split(seqname, "__"))[[1]], protcat = unlist(str_split(seqname, "__"))[[2]])
cands <- data.frame(cands, stringsAsFactors = FALSE)

#Only need those cases where a species has more than one candidate for that protein category.
cands %<>%
  group_by(species, protcat) %>%
  add_tally() %>%
  ungroup() %>%
  filter(n > 1)

#Creating a grouping variable for looping.
cands %<>% mutate(grp = paste0(species, "__", protcat))
cands %<>% arrange(desc(n))


#I will first have to write these out as their own files
#and then have msa::msa() read them in again.
fadir <- paste0(outpath, "/", "to_aln")
if(!dir.exists(fadir)){dir.create(fadir, recursive = TRUE)}

#Looping and creating the MSAs and visualizing them.
for(i in unique(cands$grp)){
  #i <- unique(cands$grp)[1]
  cur <- cands %>% filter(grp == i)
  cur %<>% transmute(seq = paste0(">", seqname, "\n", seq))
  outfile <- paste0(fadir, "/", i, ".fasta")
  write.table(cur, file = outfile, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
}
rm(cur)


##MSA GENERATION##
#MSA output directory
msadir <- paste0(outpath, "/", "msa")
if(!dir.exists(msadir)){dir.create(msadir, recursive = TRUE)}
#Listing files in fadir, and submitting these to msa:msa()
fafiles <- list.files(fadir, pattern = ".fasta", full.names = TRUE)

for(file in fafiles){
  
  #file <- fafiles[1]
  fname <- str_replace(basename(file), "\\.fasta", "")
  
  cat("Processing ", fname, "!!")
  
  alnfile <- paste0(msadir, "/", fname, "_aln.fasta")
  alnpdf <- paste0(msadir, "/", fname, "_aln.pdf")
  
  cur <- readAAStringSet(file)
  #if(!file.exists(alnfile)){
  aln <- msa::msa(cur, method = "ClustalOmega", type = "protein", 
                  order = "aligned", substitutionMatrix = "BLOSUM65")
  print(aln, show="complete")
  
  #First saving the alignment itself.
  msaPrettyPrint(aln, output = "asis", alFile = alnfile, showConsensus = "none",
                 askForOverwrite = FALSE, verbose = TRUE)
  #}
  
  
  #Visualizing quickly with seqvisr::msavisr.
  #Doing this before pretty printing the MSA because
  #that involves irreversible modifications.
  pltfile <- paste0(msadir, "/", fname, "_plt.pdf")
  #Want to put the longest file on top.
  if(!file.exists(pltfile)){
    snames <- data.frame(name = aln@unmasked@ranges@NAMES)
    snames %<>% mutate(len = str_extract(name, "(?<=len_).*$"))
    ref <- snames %>% filter(len == max(len)) %>% select(name)
    ref <- as.character(ref)
    plt <- msavisr(mymsa = alnfile, myref = ref, cbfcols = FALSE)
    plt <- plt + theme(text = element_text(size = 16))
    ggplot2::ggsave(filename = pltfile, plot = plt, device = "pdf", dpi = 1200, 
                    width = 20, height = 10)
  }
  
  if(!file.exists(alnpdf)){
    ##Now doing the pretty printing with .
    #Will put the sequence names in the legend via extra code for a caption.
    cnames <- str_replace_all(aln@unmasked@ranges@NAMES, "^.*__.*__", "")
    cnames <- str_replace_all(cnames, "_", " ")
    cnames <- paste0("S", 1:length(cnames), " - ", cnames)
    cnames <- c(str_replace_all(fname, "_+", " "), cnames)
    cnames <- paste0(cnames, collapse = "\\\\")
    capt <- paste0("\\showcaption{", cnames, "}")
    
    #Renaming for pretty printing.
    aln@unmasked@ranges@NAMES <- paste0("S", 1:length(aln@unmasked@ranges@NAMES))
    
    #MSA pretty-print does not store the files at the right location.
    #So setting the output directory for this manually.
    setwd(msadir)
    
    #Pretty printing.
    msaPrettyPrint(aln, output = "pdf", showNames = "left", file = alnpdf, consensusColors = "ColdHot",
                   showConsensus = "none", showLogo = "none", askForOverwrite = FALSE, verbose = TRUE, 
                   paperWidth = 8.3, paperHeight = 11.7, margins = c(0.75, 0.75), showNumbering = "none",
                   shadingMode = "similar", shadingModeArg = 70, furtherCode = capt)
    
    #Resetting the setwd() path to what it was before in this script.
    setwd(outpath)
  }
  
  cat("Done processing ", fname, "!!")
  
}



##Paralog vs. genetic variation estimation##
#All we can really do is just look at the
#pairwise sequence similarity distribution.
#Path where the alignments are stored.
alnfiles <- list_files(msadir, pattern = "*_aln.fasta", full.names = TRUE)
alnfiles


for(aln in alnfiles){
  
  #aln <- alnfiles[[34]]
  
  #Reading in the alignment fasta file.
  alndat <- seqvisr::fastodf(path = aln, seqtype = "AA")
  #Extracting length.
  alndat %<>% mutate(len = str_extract(seqname, "(?<=len_).*$"))
  #Extracting the original sequence.
  alndat %<>% mutate(origseq = str_replace_all(seq, "-", ""))
  
  #Pairwise identity comparisons
  #To create the pairwise combinations, throwing in
  #the sequence identifier and sequence as a single string first.
  locseqs <- alndat %>%
    mutate(origseq = paste0(seqname, "____", origseq)) %>%
    distinct(origseq)
  locseqs <- as.vector(locseqs$origseq)
  
  #Creating a data.frame of all unique pairwise combinations
  idendat <- gtools::combinations(n = length(locseqs), r = 2, v = locseqs, repeats.allowed = FALSE)
  idendat <- data.frame(idendat)
  names(idendat) <- c("query", "target")
  rm(locseqs)

  #Putting the sequence identifiers and sequences themselves
  #in separate columns.
  idendat %<>%
    rowwise() %>%
    mutate(qname = unlist(str_split(query, "____"))[[1]],
           qseq = unlist(str_split(query, "____"))[[2]],
           tname = unlist(str_split(target, "____"))[[1]],
           tseq = unlist(str_split(target, "____"))[[2]],
           .keep = "unused")
  idendat <- data.frame(idendat, stringsAsFactors = FALSE)
  
  #Pairwise identities.
  idendat %<>%
    rowwise() %>%
    mutate(pid = Biostrings::pid(Biostrings::pairwiseAlignment(qseq, tseq, substitutionMatrix = data(BLOSUM62))), 
           conseq = consensusString(consensusMatrix(pairwiseAlignment(qseq, tseq, substitutionMatrix = data(BLOSUM62)))))
  
  if(aln == alnfiles[[1]]){
    iddat <- idendat
  } else{
    iddat <- bind_rows(iddat, idendat)
  }
  
}
rm(cur, idendat, alndat)


#Plotting distribution of pairwise sequence identities.
#Binning based on this https://stackoverflow.com/a/43747717
plt <- ggplot(data = iddat, mapping = aes(x = pid)) + 
  geom_histogram(binwidth = 2, boundary = 2, fill = "darkred", colour = "black") + 
  scale_x_continuous(breaks = seq(0, 100, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) + 
  labs(y = "Count of pairwise sequence comparisons", x = "Percentage sequence identity") + 
  theme_classic() +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 20, color = "black"))
pltout <- paste0(outpath, "/", "rscript_ccseqcomp_pid_distribution.pdf")
ggsave(filename = pltout, plot = plt, dpi = 600, height = 10, width = 20)



