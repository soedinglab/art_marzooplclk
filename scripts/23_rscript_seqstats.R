#Paper's main figures and tables for the de novo assembly and annotation steps.
#Under mypath (see below), needs the general_stats/ directory
#and rscript_annotscmp/ dir with the requisite statistics files.
#These statistics files can be obtained by running the appropriate scripts
#in the scripts section, and then copying them over to the aforementioned
#directories.


rm(list = ls())

#For tables.
library(knitr)
library(kableExtra)

#For plotting.
library(showtext)
library(grDevices)
library(ggplot2)
library(ggpmisc) #To do sub-plot tables and so forth.

#General-purpose.
library(data.table)
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)


#Main path.
mypath <- "/path/to/outputs" #SET THIS PATH HERE.
setwd(mypath)

outdir <- paste0(mypath, "/", "rscript_figstabs")
if(!dir.exists(outdir)) { dir.create(outdir) }

#----------------------------------------------------------------

##STATISTICS FOR RNA SEQUENCING AND READ FILTERING##


#Reading in general stats data.
gsdir <- paste0(mypath, "/", "general_stats")
#Checking the files in there.
list.files(gsdir)


#Bowtie2 stats.
bt2dat <- fread(list.files(gsdir, pattern = "^bowtie2", full.names = TRUE), header = FALSE)

#Read counts.
readstab <- lapply(list.files(gsdir, pattern = "^seqkit", full.names = TRUE), function(x){
  df <- fread(x, header = TRUE)
  df$filename <- x
  return(df)
})
readstab <- bind_rows(readstab)

#BUSCO stats.
busdat <- fread(list.files(gsdir, pattern = "^busco", full.names = TRUE), header = FALSE, sep = "\t")



#Reading in the annotation statistics data.
annotsdir <- paste0(mypath, "/", "rscript_annotscmp")
list.files(annotsdir)
annotsdat <- fread(list.files(annotsdir, pattern = "_stats.tsv$", full.names = TRUE), header = TRUE)


#----------------------------------------------------------------



#Will start by munging the reads data table.

#Retaining only necessary columns.
readstab %<>% select(-c(format, type))

#Extracting the processing step data (e.g., raw, or after fastp)
readstab %<>% mutate(filename = basename(filename))
readstab %<>% 
  mutate(procstep = str_extract(basename(filename), "(?<=seqkit_).*(?=_stats)")) %>%
  select(-filename)

#Extracting the paired end read from the file column, and setting up the sample name.
readstab %<>% 
  mutate(read = str_extract(file, "(?<=R)[1-2]"), 
         sample = str_extract(file, "^[A-Za-z]+_[a-z]+")) %>%
  select(-file)

#The values have been presented for each read mate separately.
#Averaging this out.
#sum_len's summary should be SUM of the two because it's just
#the tota number of nucleotides.
readstab %<>% 
  select(-read) %>%
  group_by(sample, procstep) %>%
  summarise(num_seqs = as.numeric(mean(num_seqs)), 
            sum_len = as.numeric(sum(sum_len)), 
            min_len = as.numeric(mean(min_len)), 
            avg_len = as.numeric(mean(avg_len)), 
            max_len = as.numeric(mean(max_len)))
readstab <- data.frame(readstab, stringsAsFactors = FLASE)

#Write this table out as it is. It is a supplement.
outfile <- "all_read_processing_stats.csv"
outfile <- paste0(outdir, "/", outfile)
write.table(readstab, file = outfile, sep = ",", quote = TRUE, row.names = FALSE)
rm(outfile)

#Continue processing for the LATEX table.
#I will only show number of seqs before and after all processing.
readstab %<>% select(-c(min_len, avg_len, max_len, sum_len))

#I need to order the procstep properly.
readstab$procstep <- factor(readstab$procstep, levels = c("raw", "rcor", "fp", "kra2", "smr"))

#I only want to show the reads counts and what not as before processing and after processing.
#So only need rows corresponding to raw and smr.
readstab %<>% filter(procstep %in% c("raw", "smr"))




#----------------------------------------------------------------


#READS COUNTS PLOT FOR PRESENTATION.



#Plot version
readsplt <- readstab

xloc_init <- readsplt %>% filter(procstep == "raw") %>% summarize(mean = round(mean(num_seqs)/1000000))
xloc_init <- as.numeric(xloc_init)
xloc_fin <- readsplt %>% filter(procstep == "smr") %>% summarize(mean = round(mean(num_seqs)/1000000))
xloc_fin <- as.numeric(xloc_fin)


plt <- readsplt %>% 
  arrange(sample) %>%
  mutate(sample = str_replace(sample, "grp$", ""), 
         sample = str_replace(sample, "_", " "), 
         sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample)) %>%
  ggplot(mapping = aes(x = num_seqs/1000000, y = reorder(sample, desc(sample)), fill = procstep)) + 
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) + 
  scale_x_continuous(labels = scales::comma, breaks = scales::breaks_extended(n = 10)) + 
  #VLINE+TEXT FOR INITIAL READS.
  geom_vline(aes(xintercept = xloc_init), linetype = "dashed", color = "#ff3f0d", size = 1) + 
  #This subset below in geom_text() is taking only the first row, 
  #or that 50000 label would get plotted nrow() number of times.
  geom_text(data = readsplt[1,],
            aes(x = xloc_init, label = paste0("Mean = ", xloc_init, " M"), y = 9, vjust = -1.0), 
            angle = 270, size = 6, fontface = "bold", 
            color = "black") + 
  #VLINE+TEXT FOR FINAL READS.
  geom_vline(aes(xintercept = xloc_fin), linetype = "dashed", color = "#009b57", size = 1) + 
  geom_text(data = readsplt[1,],
            aes(x = xloc_fin, label = paste0("Mean = ", xloc_fin, " M"), y = 9, vjust = -1.0), 
            angle = 270, size = 6, fontface = "bold", 
            color = "black") + 
  theme_classic() + 
  xlab("Million reads") + 
  ylab("") + 
  scale_fill_manual(values = c("#ff3f0d", "#009b57"),
                    name = "", 
                    labels = c("Before processing", "After processing")) + 
  theme(legend.position = "top", 
        text = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(face = "bold.italic", size = 20)) + 
  theme(plot.title = element_text(hjust = 0.5)) #+ 
  #labs(title = "Read counts")

plt

#Saving PNG.
#outfile <- "all_reads_counts_plot.png"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "png")
#Saving PDF.
outfile <- "all_reads_counts_plot.pdf"
outfile <- paste0(outdir, "/", outfile)
ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "pdf")
#Saving SVG.
#outfile <- "all_reads_counts_plot.svg"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "svg")


rm(outfile, readsplt, xloc_fin, xloc_init)


#----------------------------------------------------------------

#READS COUNTS + ASSEMBLY STATISTICS TABLE.

#Continuing to process readstab.
#I want to pivot those columns wider and rename them.
readstab %<>% pivot_wider(names_from = procstep, values_from = num_seqs)

#Need to drop the grp from the samp column values.
readstab %<>% mutate(sample = str_replace(sample, "grp$", ""))

#I want to add in the transcriptome assembly stats from annotsdat here.
#That will be presented as its own table.
assemtab <- annotsdat %>% select(c(sample, nseqs_initial, nseqs_final, n50_initial, n50))

#Merging
assemtab <- merge(readstab, assemtab, by = "sample", all.x = TRUE, all.y = TRUE)

#Prepping some columns
#Mainly adding percentage representations to columns that are "derived"
#from other columns by filtering (e.g., initial and final assembly)
assemtab %<>%
  mutate(sample = str_replace(sample, "_", " "), 
         sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample)) %>%
  mutate(smr = paste0(smr, " (", round(smr/raw*100, 1), ")")) %>%
  mutate(nseqs_final = paste0(nseqs_final, " (", round(nseqs_final/nseqs_initial*100, 1), ")"))
  
#Writing out this latex table.
genft <- "Numbers in parentheses indicate percentage value of that column w.r.t. its counterpart."
assemtab_out <- assemtab %>%
  rename("Organism" = sample, "Raw" = raw, "Processed" = smr, 
         "Initial" = nseqs_initial, "Final" = nseqs_final,
         "N50 (initial)" = n50_initial, "N50 (final)" = n50) %>%
  kbl(., format = "latex", booktabs = T, align = "lcccccc", linesep = "") %>%
  kable_styling(position = "center", latex_options = c("scale_down")) %>%
  row_spec(row = 0, italic = FALSE, bold = TRUE) %>%
  column_spec(column = 1, italic = TRUE) %>%
  add_header_above(c(" " = 1, "Number of reads" = 2, "Transcripts in assembly" = 2, "N50" = 2), bold = TRUE) %>%
  footnote(general = genft, threeparttable = T)

assemtab_out

#Using save_kable with keep_tex = TRUE as a proxy to write the table out.
#Does not actually create the PDF since the GWDG R installation is missing
#dependencies.
#But the tex file does get written. This is all that's needed.
cat(assemtab_out, file = paste0(outdir, "/", "all_assem_stats_tab.tex"))
save_kable(assemtab_out, 
           file =  paste0(outdir, "/", "all_assem_stats_tab.pdf"), 
           keep_tex = TRUE)
#save_kable throws an error, this one to be specific https://github.com/haozhu233/kableExtra/issues/118
#None of those solutions work. But that's okay, I just need that PDF, I have the tex file anyway.

rm(assemtab_out, genft)


#----------------------------------------------------------------

#ASSEMBLY STATISTICS FIGURE FOR PRESENTATION
assemtab <- annotsdat %>% select(c(sample, nseqs_initial, nseqs_final, n50_initial, n50))

#Pivoting longer..
#Need a spec for this.
myspec <- tribble(
  ~.name, ~.value, ~seq_set,
  "nseqs_initial", "assem", "initial",
  "nseqs_final", "assem", "final",
  "n50_initial", "n50", "initial",
  "n50", "n50", "final"
)
assemtab %<>% pivot_longer_spec(spec = myspec)
rm(myspec)

#Just trying out a "long" formatted table.
#Doesn't look comfortable in latex, although the table itself is nice.
# assemtab %>%
#   kbl(., format = "latex", booktabs = T, align = "lccc", linesep = "") %>%
#   kable_styling(position = "center", latex_options = c("scale_down")) %>%
#   row_spec(row = 0, italic = FALSE, bold = TRUE) %>%
#   column_spec(column = 1, italic = TRUE) %>%
#   collapse_rows(columns = 1, valign = "top", latex_hline = "none", row_group_label_position = "identity")

#Plot.

#Need to pivot the table longer once again to get everything into the column format for ggplot.
assemtab %<>% pivot_longer(cols = c(assem, n50), names_to = "stat", values_to = "count")

#Setting levels.
#assemtab %<>% mutate(seq_set = factor(seq_set, levels = c("initial", "final")))

#Precalculating datasets for geom_vline()s and corresponding geom_texts() that will
#indicate the mean number of transcripts in the initial and final assemblies
#in the facet for that in the plot below.
vdfas <- assemtab %>% filter(stat == "assem")
xloc_vdfas_init <- vdfas %>% filter(seq_set == "initial") %>% summarize(mean = round(mean(count)))
xloc_vdfas_init <- as.numeric(xloc_vdfas_init)
xloc_vdfas_fin <- vdfas %>% filter(seq_set == "final") %>% summarize(mean = round(mean(count)))
xloc_vdfas_fin <- as.numeric(xloc_vdfas_fin)

#vline + text for n50.
vdfnf <- assemtab %>% filter(stat == "n50")
xloc_vdfnf_init <- vdfnf %>% filter(seq_set == "initial") %>% summarize(mean = round(mean(count)))
xloc_vdfnf_init <- as.numeric(xloc_vdfnf_init)
xloc_vdfnf_fin <- vdfnf %>% filter(seq_set == "final") %>% summarize(mean = round(mean(count)))
xloc_vdfnf_fin <- as.numeric(xloc_vdfnf_fin)


#PLOTTING
plt <- assemtab %>% 
  mutate(sample = str_replace(sample, "_", " "), 
         sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample),
         sample = reorder(sample, desc(sample))) %>%
  group_by(sample, stat) %>%
  #arrange(ifelse(stat == "assem", rev(desc(seq_set)), seq_set), .by_group = TRUE) %>%
  #Arranging by desc(count) ensures that the value in the overlapping bars with the
  #lower value gets plotted to the front.
  arrange(desc(count), .by_group = TRUE) %>%
  #arrange(seq_set, .by_group = TRUE) %>%
  #ungroup() %>%
  ggplot(aes(x = count, y = sample, fill = seq_set)) + 
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) + 
  scale_x_continuous(labels = scales::comma, breaks = scales::breaks_extended(n = 10)) + 
  scale_fill_manual(values = c("#009b57", "#ff3f0d"),
                    name = "Assembly", , 
                    labels = c("Final", "Initial")) + 
  facet_wrap(~stat, scales = "free_x", 
             strip.position = "bottom", 
             labeller = as_labeller(c(assem = "Number of transcripts", n50 = "N50 sequence length (nt)") )) + 
  #VLINE+TEXT FOR INITIAL assembly mean number.
  geom_vline(data = vdfas, aes(xintercept = xloc_vdfas_init), linetype = "dashed", color = "#ff3f0d", size = 1) + 
  #This subset below in geom_text() is taking only the first row, 
  #or that 50000 label would get plotted nrow() number of times.
  geom_text(data = vdfas[1,],
            aes(x = xloc_vdfas_init, label = paste0("Mean = ", xloc_vdfas_init), y = 9, vjust = -1.0), 
            angle = 270, size = 6, fontface = "bold", 
            color = "black") + 
  #VLINE+TEXT FOR FINAL assembly mean number.
  geom_vline(data = vdfas, aes(xintercept = xloc_vdfas_fin), linetype = "dashed", color = "#009b57", size = 1) + 
  geom_text(data = vdfas[1,],
            aes(x = xloc_vdfas_fin, label = paste0("Mean = ", xloc_vdfas_fin), y = 9, vjust = -1.0), 
            angle = 270, size = 6, fontface = "bold", 
            color = "black") + 
  #VLINE+TEXT FOR N50 INITIAL.
  geom_vline(data = vdfnf, aes(xintercept = xloc_vdfnf_init), linetype = "dashed", color = "#ff3f0d", size = 1) + 
  geom_text(data = vdfnf[1,],
            aes(x = xloc_vdfnf_init, label = paste0("Mean = ", xloc_vdfnf_init), y = 9, vjust = -1.0), 
            angle = 270, size = 6, fontface = "bold", 
            color = "black") + 
  #VLINE+TEXT FOR N50 FINAL.
  geom_vline(data = vdfnf, aes(xintercept = xloc_vdfnf_fin), linetype = "dashed", color = "#009b57", size = 1) + 
  geom_text(data = vdfnf[1,],
            aes(x = xloc_vdfnf_fin, label = paste0("Mean = ", xloc_vdfnf_fin), y = 9, vjust = -1.0), 
            angle = 270, size = 6, fontface = "bold", 
            color = "black") + 
  theme_classic() + 
  theme(strip.placement = "outside", 
        strip.background = element_blank(),
        legend.position = "top",
        text = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(face = "bold.italic", size = 20), 
        axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(hjust = 0.5), legend.justification = "center") +
  guides(fill = guide_legend(reverse = TRUE)) + 
  #theme_classic() + 
  ylab("") + 
  xlab("") #+
  #labs(title = "Assembly comparison")

plt

#Saving PNG.
#outfile <- "all_assem_stats_plot.png"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "png")
#Saving PDF.
outfile <- "all_assem_stats_plot.pdf"
outfile <- paste0(outdir, "/", outfile)
ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "pdf")
#Saving SVG.
#outfile <- "all_assem_stats_plot.svg"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "svg")


rm(outfile, plt, vdfas, vdfnf, xloc_vdfas_fin, xloc_vdfas_init, xloc_vdfnf_fin, xloc_vdfnf_init)

rm(assemtab)




#----------------------------------------------------------------


#ANNOTATIONS FIGURE FOR PRESENTATION.

#For this, I want nseqs_final, mean_protlen, swissprot_annots, eggnogmapper_annots, miss_annots.
annotstab <- annotsdat %>% select(c(sample, mean_protlen, nseqs_final, swissprot_annots, eggnogmapper_annots, miss_annots))

#Setting up the sample names.
annotstab %<>% mutate(sample = str_replace(sample, "_", " "), 
                     sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample),
                     sample = reorder(sample, desc(sample)))

#Pivoting longer.
annotstab %<>% pivot_longer(cols = -sample, names_to = "stat", values_to = "count")

#Setting order for the stat column.
mylvls <- c("miss_annots", "swissprot_annots", "eggnogmapper_annots", "nseqs_final", "mean_protlen")
#mylvls <- c("mean_protlen", "nseqs_final", "swissprot_annots", "eggnogmapper_annots", "miss_annots")
annotstab %<>% mutate(stat = factor(stat, levels = mylvls))



#I will show the mean values for annotations as a table within the plot.
meantab <- annotstab %>% 
  filter(stat != "mean_protlen") %>%
  group_by(stat) %>% 
  summarize(mean = round(mean(count))) %>%
  mutate(meanperc = round(mean/mean[4]*100, 0)) %>%
  mutate(mean = ifelse(meanperc == 100, paste0(mean, " (", meanperc, "%)"), paste0(mean, "  (", meanperc, "%)"))) %>%
  select(-meanperc) %>%
  dplyr::arrange(-dplyr::row_number()) %>%
  mutate(stat = case_when(
    stat == "nseqs_final" ~ "Total sequences", 
    stat == "eggnogmapper_annots" ~ "eggnog-mapper annotations", 
    stat == "swissprot_annots" ~ "Swiss-Prot annotations",
    stat == "miss_annots" ~ "Unannotated")) %>%
  #Padding stat column values if not they get center aligned in inset table
  mutate(stat = str_pad(stat, width = max(nchar(stat)), side = "right"), 
         mean = str_pad(mean, width = max(nchar(mean)), side = "left")) %>%
  rename(" " = stat, "Mean counts" = mean)

meantab<- tibble(x = max(annotstab$count)+4000, 
                 y = round(length(unique(annotstab$sample))+0.5, 0),
                 #y = 0, 
                 tb = list(meantab))


#Plotting.
plt <- annotstab %>%
  group_by(sample) %>%
  arrange(stat) %>%
  ggplot(mapping = aes(x = count, y = sample, fill = stat)) +
  geom_bar(position = "identity", stat = "identity", alpha = 0.5, 
           data = annotstab[annotstab$stat == "nseqs_final", ]) + 
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8, 
           data = annotstab[!(annotstab$stat %in% c("nseqs_final", "mean_protlen")),]) + 
  geom_text(data = annotstab[annotstab$stat == "mean_protlen", ], 
            mapping = aes(y = sample, 
                          x = annotstab[annotstab$stat == "nseqs_final", ]$count, 
                          label = round(count, 0)), 
            hjust = 1.2, fontface = "bold", size = 10) + 
  scale_fill_manual(breaks = c("eggnogmapper_annots", "swissprot_annots", "miss_annots", "nseqs_final"), 
                    values = c("#33a02c", "#1f73b2", "#ff3f0d", "gray", "black"),
                    name = " ", 
                    labels = c("eggnog-mapper annotations", "Swiss-Prot annotations", "Unannotated", "Total sequences")) + 
  scale_x_continuous(labels = scales::comma, breaks = scales::breaks_extended(n = 10)) +
  theme_classic() + 
  theme(legend.position = "top",
        legend.text = element_text(size = 16),
        text = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(face = "bold.italic", size = 20), 
        #axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(hjust = 0.5), legend.justification = "center") +
  #guides(fill = guide_legend(nrow = 2)) + 
  #theme_classic() + 
  ylab("") + 
  xlab("") +
  geom_table(data = meantab, 
             aes(x = x, y = y, label = tb),
             table.theme = ttheme_gtminimal(base_family = "mono", base_size = 20)) +
  annotate(geom = "text", parse = TRUE,
           label = 'bold("Inset numbers in gray bars:\nMean protein lengths")', 
           x = max(annotstab$count)-4000, y = 1, size = 8) #+ 
  #labs(title = "Functional annotation metrics")

plt

#Saving PNG.
#outfile <- "all_annots_stats_plot.png"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "png")
#Saving PDF.
outfile <- "all_annots_stats_plot.pdf"
outfile <- paste0(outdir, "/", outfile)
ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = cairo_pdf)
#Saving SVG.
#outfile <- "all_annots_stats_plot.svg"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "svg")


rm(outfile, plt, meantab, mylvls)

rm(annotstab)



#----------------------------------------------------------------


#ANNOTATIONS TABLE.

#For this, I want nseqs_final, mean_protlen, swissprot_annots, eggnogmapper_annots, miss_annots.
annotstab <- annotsdat %>% select(c(sample, mean_protlen, nseqs_final, swissprot_annots, eggnogmapper_annots, miss_annots))

#Setting up the sample names.
annotstab %<>% mutate(sample = str_replace(sample, "_", " "), 
                      sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample),
                      sample = reorder(sample, desc(sample)))


#Rounding off that mean_protlen value.
annotstab %<>% mutate(mean_protlen = round(mean_protlen, 0))

#Getting percentages and adding them to the respective columns.
annotstab %<>%
  mutate(swissprot_annots = paste0(swissprot_annots, " (", round(swissprot_annots/nseqs_final*100, 0), "%)"), 
         eggnogmapper_annots = paste0(eggnogmapper_annots, " (", round(eggnogmapper_annots/nseqs_final*100, 0), "%)"), 
         miss_annots = paste0(miss_annots, " (", round(miss_annots/nseqs_final*100, 0), "%)")  )

#Rearranding columns a bit.
annotstab %<>% select(c(sample, nseqs_final, mean_protlen, eggnogmapper_annots, swissprot_annots, miss_annots))

#Writing out this latex table.
genft <- "Percentage values calculated as fraction of total sequences from column 2."
annotstab_out <- annotstab %>%
  rename("Organism" = sample, "Num. proteins" = nseqs_final, "Mean length" = mean_protlen, 
         "eggnog-mapper" = eggnogmapper_annots, "Swiss-Prot" = swissprot_annots, "Unannotated" = miss_annots) %>%
  kbl(., format = "latex", booktabs = T, align = "lcccccc", linesep = "") %>%
  kable_styling(position = "center", latex_options = c("scale_down")) %>%
  row_spec(row = 0, italic = FALSE, bold = TRUE) %>%
  column_spec(column = 1, italic = TRUE) %>%
  add_header_above(c(" " = 3, "Number of proteins annotated" = 3), bold = TRUE) %>%
  footnote(general = genft, threeparttable = T)

annotstab_out

#Using save_kable with keep_tex = TRUE as a proxy to write the table out.
#Does not actually create the PDF since the GWDG R installation is missing
#dependencies.
#But the tex file does get written. This is all that's needed.
cat(annotstab_out, file = paste0(outdir, "/", "all_annots_stats_tab.tex"))
save_kable(annotstab_out, 
           file =  paste0(outdir, "/", "all_annots_stats_tab.pdf"), 
           keep_tex = TRUE)
#save_kable throws an error, this one to be specific https://github.com/haozhu233/kableExtra/issues/118
#None of those solutions work. But that's okay, I just need that PDF, I have the tex file anyway.

rm(annotstab_out, annotstab)


#----------------------------------------------------------------


#ASSEMBLY QUALITY METRICS.


#BUSCO TABLE FOR SUPPLEMENT.
bustab <- busdat

#Can drop that third column.
bustab %<>% select(-V3)

#Renaming columns.
names(bustab) <- c("filename", "count")

#Need to extract the BUSCO step and sample name from the filenames first.
bustab %<>% mutate(filename = str_replace(filename, "^/cbscratch/.*/busco_", ""))

#Can drop the entire short_summary.*$ string 'cause the 
#necessary info is embedded in the rest of the filename string.
bustab %<>% mutate(filename = str_replace(filename, "/short.*$", ""))

#Now I can extract the procstep (e.g., Trinity or rfilt) and the samplename + busco set strings.
bustab %<>% 
  rowwise() %>%
  mutate(procstep = unlist(str_split(filename, "/"))[1], 
         sample = unlist(str_split(filename, "/"))[2])
bustab <- data.frame(bustab, stringsAsFactors = FALSE)

bustab %<>% select(-filename)

#Extracting the busco dataset from the sample column,
#and cleaning up the sample column.
bustab %<>% 
  mutate(busset = str_extract(sample, "(?<=_busco_)[a-z]+"), 
         sample = str_extract(sample, "^.*(?=_busco)"))

#Cleaning up the count column.
bustab %<>% mutate(count = str_replace_all(count, "[\\[\\]]+", ","))
bustab %<>% mutate(count = str_replace_all(count, ",,", ","))


#So in this data I have BUSCO scores for:
unique(bustab$procstep)
#protclu_td = trinity, intermediate, prot
#rfilt = trinity, final, nuc
#rscript_lenfilt = trinity, final, prot
#spades = spades, initial, nuc
#td_trinity = trinity, initial, prot
#trinity = trinity, initial, nuc

#Put spades and trin in a separate table for later.
spadtrin <- bustab %>% filter(procstep %in% c("spades", "trinity"))

#Dropping non-intial and final assembly BUSCOS.
bustab %<>% filter(procstep %in% c("trinity", "td_trinity", "rfilt", "rscript_lenfilt"))

#Combining the busset and ngenes columns since they represent facet panels
bustab %<>% 
  mutate(busset = case_when(
    busset == "euk" ~ "BUSCO Eukaryota", 
    busset == "arth" ~ "BUSCO Arthropoda"))

#Will extract the BUSCO dataset's gene count number and attach it to busset
bustab %<>% mutate(ngenes = str_extract(count, "(?<=n\\:)[0-9]+"))
bustab %<>% mutate(busset = paste0(busset, " (n = ", ngenes, ")"))
bustab %<>% mutate(count = str_replace(count, ",n.*$", ""))

#Fixing the procstep names.
bustab %<>%
  mutate(procstep = case_when(
    procstep == "rfilt" ~ "Final assembly",
    procstep == "rscript_lenfilt" ~ "Final proteome",
    procstep == "trinity" ~ "Initial assembly",
    procstep == "td_trinity" ~ "Initial proteome"
  ))

#Fixing sample names.
bustab %<>% mutate(sample = str_replace(sample, "grp", ""),
                   sample = str_replace(sample, "_", " "), 
                   sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample),
                   sample = reorder(sample, desc(sample)))

#Rearranging and renaming columns.
bustab %<>% select(c(sample, procstep, busset, count))


#Writing out this table.
#Write this table out as it is. It is a supplement.
outfile <- "all_busco_stats.csv"
outfile <- paste0(outdir, "/", outfile)
write.table(bustab, file = outfile, sep = ",", quote = TRUE, row.names = FALSE)
rm(outfile)




#----------------------------------------------------------------


#CONTINUING BUSCO TABLE PROCESSING FOR INTEGRATION WITH BOWTIE SCORES
#AND FOR SUBSEQUENT FIGURE OF THESE TOGETHER.


#Moving each count into its own row.
bustab %<>% separate_rows(count, sep = ",")


#Creating a new column indicating what the count value is in terms of
#BUSCO nomenclature.
bustab %<>%
  mutate(buscat = case_when(
    str_detect(count, "^C") ~ "Complete",
    str_detect(count, "^S") ~ "Single-copy (Complete)",
    str_detect(count, "^D") ~ "Duplicated (Complete)",
    str_detect(count, "^M") ~ "Missing",
    str_detect(count, "^F") ~ "Fragmented",
    str_detect(count, "^n") ~ "num_buscos",
  ))

#Setting levels for the buscat columns so that Single-copy and Duplicated will
#occur next to one another.
bustab %<>% 
  mutate(buscat = factor(buscat, levels = rev(c("Complete", "Single-copy (Complete)", "Duplicated (Complete)", 
                                                "Fragmented", "Missing"))))

#Extracting the percentages from the counts column.
bustab %<>% mutate(count = as.numeric(str_extract_all(count, "[0-9\\.]+")))




#----------------------------------------------------------------


#PREPARING BOWTIE2 DATA.
#THIS WILL BE PLOTTED ALONG WITH THE BUSCO DATA.

bt2tab <- bt2dat

#Can drop those third-fifth columns.
bt2tab %<>% select(-c(V3, V4, V5))

#Renaming columns.
names(bt2tab) <- c("filename", "count")

#Need to extract the Bowtie2 step and sample name from the filenames first.
bt2tab %<>% mutate(filename = str_replace(filename, "^/cbscratch/.*/bowtie2_", ""))
bt2tab %<>% mutate(filename = str_extract(filename, "^[A-Za-z]+\\/[A-Za-z_]+"))
bt2tab %<>% 
  rowwise() %>% 
  mutate(procstep = unlist(str_split(filename, "\\/"))[1], 
         sample = unlist(str_split(filename, "\\/"))[2])

#Dropping filename and extracting the numeric in the count column.
bt2tab %<>% 
  select(-filename) %>%
  mutate(count = as.numeric(str_replace(count, "%", "")))

#Filtering to retain only procstep of interst.
bt2tab %<>% filter(procstep %in% c("trinity", "rfilt"))

#Renaming procstep in-line with bustab.
bt2tab %<>%
  mutate(procstep = case_when(
    procstep == "rfilt" ~ "Final assembly",
    procstep == "trinity" ~ "Initial assembly"
  ))

#Fixing sample names.
bt2tab %<>% mutate(sample = str_replace(sample, "grp", ""),
                   sample = str_replace(sample, "_", " "), 
                   sample = ifelse(str_detect(sample, "sp$"), paste0(sample, "."), sample),
                   sample = reorder(sample, desc(sample)))

#Rearranging columns.
bt2tab %<>% select(c(sample, procstep, count))

#Writing out this table.
#Write this table out as it is. It is a supplement.
outfile <- "all_bowtie2_stats.csv"
outfile <- paste0(outdir, "/", outfile)
write.table(bt2tab, file = outfile, sep = ",", quote = TRUE, row.names = FALSE)
rm(outfile)



#----------------------------------------------------------------


#INTEGRATING BUSCO AND BOWTIE2 TABLES FOR PLOTTING.

#Need to add a buscat and busset column each to bt2tab.
#busset will have to be set to Bowtie2 to facet properly.
#buscat will be set to "Read support".
bt2tab %<>% mutate(buscat = "Read support", busset = "Bowtie2 read support")

#Rearranging bt2tab.
bt2tab %<>% select(c(sample, procstep, busset, count, buscat))

#Rowbinding the two data.frames.
qualdat <- bind_rows(bt2tab, bustab)




#----------------------------------------------------------------





#PLOTTING QUALITY ASSESSMENT.

#Looking at the procsteps I have.
unique(qualdat$procstep)

#As the BUSCO completeness scores are roughly equivalent at both the
#nucleotide and protein levels, I will just go with the nucleotide
#level stuff (initial and final "assembly"). This also makes it
#easier to plot the bowtie2 scores alongside which only exist
#at the nucleotide level.
qualdat %<>%
  filter(procstep %in% c("Initial assembly", "Final assembly"))

#Setting up levels.
buscat_busc_lvls <- c("Complete", "Single-copy (Complete)", "Duplicated (Complete)", "Fragmented", "Missing")
buscat_lvls <- c(buscat_busc_lvls, "Read support")
qualdat %<>% mutate(procstep = factor(procstep, levels = c("Initial assembly", "Final assembly")), 
                    buscat = factor(buscat, levels = buscat_lvls), 
                    busset = factor(busset, levels = c("BUSCO Arthropoda (n = 1013)", "BUSCO Eukaryota (n = 255)", 
                                                       "Bowtie2 read support")))
rm(buscat_busc_lvls, buscat_lvls)

#I want to add geom_vlines() to describe the respective mean values for each facet.
meantab <- qualdat %>%
  filter(buscat %in% c("Complete", "Read support")) %>%
  group_by(procstep, busset) %>%
  summarize(mean = round(median(count), 2))

#I will merge this in to qualtab and supply it as a column.
meantab %<>% data.frame() %>%
  mutate(mergecol = paste0(procstep, "__", busset)) %>%
  select(-c(procstep, busset))

#Merging.
qualdat %<>% mutate(mergecol = paste0(procstep, "__", busset))
qualdat <- merge(qualdat, meantab, by = "mergecol", all.x = TRUE)
#qualdat %<>% select(-mergecol)
#rm(meantab)

#Filtering.
qualdat %<>%
  filter(buscat != "Complete")

#Ordering samples properly.
#qualdat <- data.frame(qualdat)
qualdat %<>% 
  group_by(procstep, busset, buscat) %>%
  mutate(sample = as.character(sample)) %>%
  arrange(sample, .by_group = TRUE) %>%
  mutate(sample = reorder(sample, desc(sample))) %>%
  ungroup()


#PLOTTING.
plt <- qualdat %>%
  ggplot(aes(x = count, y = sample, fill = buscat)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) + 
  facet_grid(procstep ~ busset, switch = "y", space = "free", shrink = FALSE) + 
  geom_vline(mapping = aes(xintercept = mean), linetype = "dashed", color = "black", size = 1, alpha = 1.0) + 
  geom_text(data = filter(.data = qualdat, mergecol == meantab$mergecol[1])[1, ],
            aes(x = meantab$mean[1], label = paste0("Mean = ", meantab$mean[1]), y = 9, vjust = 1.2), 
            angle = 270, size = 6, fontface = "bold", color = "black") + 
  geom_text(data = filter(.data = qualdat, mergecol == meantab$mergecol[2])[1, ],
            aes(x = meantab$mean[2], label = paste0("Mean = ", meantab$mean[2]), y = 9, vjust = 1.2), 
            angle = 270, size = 6, fontface = "bold", color = "black") + 
  geom_text(data = filter(.data = qualdat, mergecol == meantab$mergecol[3])[1, ],
            aes(x = meantab$mean[3], label = paste0("Mean = ", meantab$mean[3]), y = 9, vjust = 1.2), 
            angle = 270, size = 6, fontface = "bold", color = "black") + 
  geom_text(data = filter(.data = qualdat, mergecol == meantab$mergecol[4])[1, ],
            aes(x = meantab$mean[4], label = paste0("Mean = ", meantab$mean[4]), y = 9, vjust = 1.2), 
            angle = 270, size = 6, fontface = "bold", color = "black") + 
  geom_text(data = filter(.data = qualdat, mergecol == meantab$mergecol[5])[1, ],
            aes(x = meantab$mean[5], label = paste0("Mean = ", meantab$mean[5]), y = 9, vjust = 1.2), 
            angle = 270, size = 6, fontface = "bold", color = "black") + 
  geom_text(data = filter(.data = qualdat, mergecol == meantab$mergecol[6])[1, ],
            aes(x = meantab$mean[6], label = paste0("Mean = ", meantab$mean[6]), y = 9, vjust = 1.2), 
            angle = 270, size = 6, fontface = "bold", color = "black") + 
  theme_classic() + 
  theme(legend.position = "top",
        legend.text = element_text(size = 20),
        text = element_text(size = 20, face = "bold"), 
        axis.text.y = element_text(face = "bold.italic", size = 14), 
        #axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(hjust = 0.5), legend.justification = "center", 
        strip.background.x = element_blank(), 
        strip.text = element_text(size = 20, face = "bold"), 
        strip.placement = "outside") +
  guides(fill = guide_legend(nrow = 1, title = "")) + 
  scale_fill_manual(name = "BUSCO gene", values = c("#33a02c", "#b2df8a", "#fc8d62", "#377eb8", "#de2d26")) +
  ylab("") + 
  xlab("Percentage")
  
plt



#Saving PNG.
#outfile <- "all_qual_stats_plot.png"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "png")
#Saving PDF.
outfile <- "all_qual_stats_plot.pdf"
outfile <- paste0(outdir, "/", outfile)
ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = cairo_pdf)
#Saving SVG.
#outfile <- "all_qual_stats_plot.svg"
#outfile <- paste0(outdir, "/", outfile)
#ggsave(outfile, plot = plt, width = 20, height = 10, dpi = 600, device = "svg")


rm(outfile, plt, meantab)

rm(qualdat)



#----------------------------------------------------------------







