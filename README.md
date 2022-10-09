RNA sequencing reveals circadian clocks may be widespread in marine zooplankton
================

- This repository contains **data** and **scripts** from the publication "RNA sequencing reveals circadian clocks may be ubiquitous in marine zooplankton".

- **Brief summary of the project:** we explored the _de novo_ assembled transcriptomes of a diverse set of marine zooplankton species to reveal the presence of transcripts coding for proteins that may together constitute functional circadian clocks. This was achieved by identifying candidate sequences via orthology to a select set of "bait" circadian clock protein sequences from a set of model organisms.

- Note: some data from this study (sequencing reads, for example) are not hosted here, and can instead be found on the NCBI under the accession code `PRJNA824716`.

- Citations:
  - Manuscript: in submission.
  - This repository: [![DOI](https://zenodo.org/badge/496647002.svg)](https://zenodo.org/badge/latestdoi/496647002)

- **Please use our issue tracker or email us for further assistance.**


## Data

- `annotations` contains compressed flat files (`*.tsv.xz`) that tabulate sequence annotations for the _de novo_ assembled transcripts. Annotations include homologs from `Swiss-Prot`, and orthologs, domains, and gene ontology terms from `EggNOG`.

- `transcriptomes` and `proteomes` contain the _de novo_ assembled transcriptome assemblies and their _in silico_ translated protein sequence sets (proteomes) respectively.

- `circadian_clock_candidates` contains `FASTA` files of sequences identified as circadian clock proteins. The `fas_by_samp` sub-directory presents these grouped by host species, and the `fas_by_type` sub-directory presents the very same sequences, but grouped by clock protein type. The `domain_structure_visualizations` sub-directory contains visualizations of the functional domains in the protein sequences, grouped together by circadian clock protein type. Finally, the `circadian_clock_candidates` directory also contains a file tabulating annotations for the candidate sequences (`cc_cand_sel_main_table.csv`). Within this table, each row corresponds to a particular candidate sequence. The columns represent different kinds of annotations and sequence-related information. The contents of the columns are as follows:

```
orthogroup - EggNOG orthogroup the sequence was assigned to. 

refseq	- sequence identifer of the bait that was used to find the candidate sequence.
matchseq - sequence identifier of the candidate.

spref - species to which the bait belongs.
spmatch	- species to which the candidate belongs.
refmatcat	- as the search strategy also allows for the bait to identify other baits that are orthologs, the matches in the table may be classified into two categories - bait vs. bait (SOIREF__SOIREF) and bait vs. candidate (SOIREF__TRINITY).

pacc_match - PANTHER DB annotation's accession for the candidate.
pdesc_match - PANTHER DB annotation's description for the candidate.
pacc_ref - PANTHER DB annotation's accession for the bait.
pdesc_ref	- PANTHER DB annotation's description for the bait.

pfacc_match - Pfam DB annotation's accession for the candidate.
pfdesc_match - Pfam DB annotation's description for the candidate.
pfloc_match - Pfam domains' locations on on the sequence for the candidate.
pfacc_ref - Pfam DB annotation's accession for the bait
pfdesc_ref - Pfam DB annotation's description for the bait.
pfloc_ref - Pfam domains' locations on on the sequence for the bait.

protcat - which type of circadian clock protein the candidate is.
mconf - binary indicator of whether the PANTHER DB sub-family annotations of the bait and the candidate are identical or non-identical (1 - identical, 0 - non-identical).

upid - UniProt identifier of the best match homolog for the candidate.
swissprot_subject - UniProt/Swiss-Prot sequence header string of the best match homolog for the candidate.
swissprot_percid - pairwise percentage identity between the sequences of the best match homolog from Swiss-Prot and the candidate.
swissprot_evalue - likelihood of this sequence from Swiss-Prot having been assigned as the best match to the candidate purely by chance (E-value; smaller the value, more unlikely that the assignment was by chance).

fas_match - protein sequence of the candidate.
fas_ref - protein sequence of the bait.

pid - pairwise percentage identity between the sequences of the bait and the candidate.

candname - sequence header of the candidate.
```


## Scripts

- The scripts used to execute the tools and perform all analyses for the paper are contained in this directory. Please refer to the workflow figure included in the "Methods" section of the manuscript for a diagrammatic overview.

- The last few R scripts together with the scripts in the `of_ref_prep` and `species_tree` also produce all plots and figures in the manuscript.

- A vast majority of the scripts here are `BASH` scripts that are used to submit jobs on a cluster.

- All scripts were executed on a high performance compute cluster running `Scientific Linux release 7.9 (Nitrogen)` with `SLURM` as the scheduler.

- For easy replication of the scripts here, it is best to organize everything within a directory structure like so (with the names indicated below):

```
newrun/
 |-jobfiles/
 |-cc_orthofinder_dbs/
 |-reads/
 |-outputs/
```

- All scripts discussed below should be placed in the `jobfiles/` directory (indicated in the directory tree above); these are some exceptions, and these will be indicated as they arise. When executed, they will place outputs in the `newrun/outputs` directory, and will place `stdout` and `stderr` data + EXECUTION scripts (see next point below) under `jobfiles/<appropriately_named_directory>`. There are a few `R` scripts in this workflow. These do not save `stderr` and `stdout` by default, but rather only print to the terminal instead. The user can optionally redirect the output to a file to create a log of the run.

- Recreating the workflow is simple once all the data and scripts and so forth are placed wherever they should be. The user must just execute the scripts in the order indicated by the numerical prefixes to the script names. I.e., `1_*_mainjob.bash` must be executed first, then `2_*_mainjob.bash`, and so forth. If there are scripts with additional prefixes (e.g., `1a_` and `1b_`), the script with the prefix occuring first in the alphabetical order must be executed first (i.e., `1a_` before `1b_`). The script names should be fairly self-explanatory. In pretty much all cases, setting the main path (via the `MYPATH` variable or its equivalent, e.g., `mypath`) to point to wherever the `newrun/` directory is should suffice.

- As these scripts were designed to run on compute nodes on a `SLURM`-based cluster, the thread counts and memory allocations are large in many cases (typically 16 cores and 124GB of RAM). The main `OrthoFinder` run uses a full complement of 256 or 128 cores + all available memory. The user will have to change these parameters accordingly for their use case. Additionally, some of the scripts by contain copy/move commands to move inputs and outputs onto--or off of--local storage on the compute nodes. These must be adjusted according to local requirements.

- The tools used are, in most cases, installed via `conda`; conda activation and deactivation invocations are therefore to be found inscribed within the job scripts. Tool version numbers are noted in the manuscript. Please modify the scripts accordingly before use should `conda` not be available.

- Most scripts here are `BASH` scripts meant to be executed on a `SLURM`-based high performance compute cluster. In most cases, they are also SUBMITTER scripts: scripts carrying a `for loop` within that cycles through the inputs at a given path, and creates EXECUTION `BASH` scripts for each input(s). The SUBMITTER also submits these EXECUTION scripts to the `SLURM` scheduler via `sbatch`. That is, most scripts will have the following generic structure:

```
#sbatch <slurm_stuff>

#This is the SUBMITTER script.
PATH="/path/to/inputs";

#Constructing the EXECUTION script using a for loop holding a template.
for SAMPLE in PATH;
do
    #Set up some stuff.
    EXECSCRIPT="/path/to/execscript/$(basename ${SAMPLE})";
    echo "
#sbatch <slurm_stuff>
#This is the EXECUTION script's template body.

#Executing tool.
mytool -i ${SAMPLE} -o /path/to/output/$(basename ${SAMPLE});
" > ${EXECSCRIPT};
    
    #Executing the exec script by submitting it to sbatch.
    sbatch ${EXECSCRIPT};

done;
```

- Non-`SLURM` scripts also exist in this workflow. None of this code is parellelized (i.e., they are not submitter scripts). These are:

```
10b_sed_td_postprocess_mainjob.bash - pass the path to the TransDecoder output directory as input, execute directly at the command line.
13_rscript_lenfilt_new_mainjob.R - execute via Rscript directly from command line, set input paths within script properly before execution.
19_rscript_of_getsois_mainjob.R - execute via Rscript, inputs must be set within the script, input is basically the path to the OrthoFinder outputs directory. This one is a bit strange in that it places its outputs under the outputs/orthofinder/ directory.
20c_awk_sed_interproscan_postprocessing_mainjob.bash - pass the path to the InterProScan outputs of 19_rscript_of_getsois_mainjob.R, execute directly from the command line via Rscript.
21a_rscript_ccsel_main.R - execute via Rscript, inputs must be set within the script (pay heed to the comments), must have 21b_rscript_ccsel_auxfunc.R alongside in the same directory to function properly.
17_rscript_annotscmp_mainjob.R - execute via Rscript, inputs must be set within the script.
23_rscript_seqstats.R and 24_rscript_ccvissum_mainjob.R - execute via Rscript, see comments in script. These two generate the final tables, plots, and figures.
```

- As can be seen from some of the script names above (e.g., `22_rscript_ccsel_main.R`), quite a few scripts use the `R` programming language (`v3.6.0` or greater) and concomitant packages. An exhaustive list of packages used in these scripts (and package versions) can be found in the file `r_packages_used.csv` in the `scripts/` directory. Accessing the scripts through `RStudio` should highlight missing packages automatically and present an option to install them automatically; this does not work for the `Bioconductor` packages nor packages from `GitHub` (e.g., `seqvisr`).

- Very important: some scripts contain additional post-processing commands built into them. For instance, the script `06_trinity_mainjob.bash`, which runs the `Trinity` de novo transcriptome assembler, has a `SeqKit` call after the main tool's execution call to filter out all assembled contigs that are shorter than 200 nucleotides (this seems to be a Trinity bug at the moment that 200 nucleotides is set as the minimum length, but shorter sequences are emitted anyway). Likewise the script that executes the `TransDecoder` sequence translation tool has quite a few post-processing commands built in also. We have not provided an overview of what post-processing is in-built into what script, and urge the user to take note of such instances by reading the code.

- Databases: most databases used have versioning bound to their tools (e.g., `SortMeRNA` and its RNA database, `Kraken2` and the `PlusPFP` DB). The only real "external" databases used are `UniProt/SwissProt` `v2021_03` (for transcriptome annotation via `MMseqs2`), and the three reference proteomes from `UniProt` used by `OrthoFinder` to find the circadian clock proteins. These are `Danaus_plexippus_UP000007151_278856.fasta`, `Drosophila_melanogaster_UP000000803_7227.fasta`, and `Mus_musculus_UP000000589_10090.fasta` (all downloaded on or before March 2021). There are copies of these in the `of_ref_prep/proteomes` directory, and do not have to be generated afresh to replicate this workflow.

- In the interest of easing the burden of reproducing the results here, the most important database resources have been provided under the directory `databases_used/` in this repository. Place the `cc_orthofinder_dbs/` directory from here directly under `newrun/` as indicated in the directory tree earlier above. The `UniProt/SwissProt` `FASTA` file is also provided here under `databases_used/` as an `xz`-compressed file. This must be decompressed (decompress with `xz -vd <filename>`) and placed at a convenient location (e.g., directly under `newrun/`), and this path then edited into line 13 of `16b_mmseqs2_annot_mainjob.bash` prior to execution.

- The scripts and data in the `scripts/of_re_prep` sub-directory in this collection of scripts need not necessarily be run. Their outputs are already present (the `of_ref_prep/final_ref_orthofinder` directory). These constitute the set of reference transcriptomes supplied to `OrthoFinder` to identify circadian clock protein candidates in the other input proteomes by means of pairwise orthology to select bait sequences (appropriately highlighted by annotations) in a select set of reference proteomes.

- The script (and data) in the `scripts/ncbi_assem_counts` directory are not a part of the analytical workflow. These were used to estimate how many of the sequenced species had genomes/transcriptomes available on NCBI.

- Similarly, the script and data in the `scripts/species_tree` directory are also not a part of the analytical workflow. These were used to generated the species tree figure included in the manuscript.


## Miscellaneous

- `supplementary_files` contains a copy of the supplementary materials submitted along with the manuscript.

- `ncbi_submission` contains some material concerning the NCBI submission for this project. Data uploaded to NCBI can be found under the accession code `PRJNA824716`.
