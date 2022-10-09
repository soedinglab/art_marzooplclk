rm(list = ls())


library(Rcpp)
library(seqinr)
library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)


#Function to return a data.frame given the path to a fasta file.
#The data.frame contains the sequence's name, sequence, and the filename (in that order).
#This function needs the seqinr package.
fastodf <- function(x){
  
  require(seqinr)
  
  curdat <- seqinr::read.fasta(file = x, as.string = TRUE, whole.header = TRUE, seqtype = "AA")
  
  curdat <- data.frame(curname = names(unlist(curdat)), curseq = unlist(curdat), stringsAsFactors = FALSE)
  curdat$filename <- rep_len(basename(x), length.out = nrow(curdat))
  
  row.names(curdat) <- NULL
  
  return(curdat)
  
}

#From here: https://stackoverflow.com/a/64404233
#String compare and extract substring function.
cppFunction('
String largest_common_substring(String str1, String str2) 
{ 
    std::string S = str1;
    std::string T = str2;
    int r = S.length();
    int n = T.length();
    std::vector<std::vector<int> > L(r , std::vector<int>(n));
    int z = 0;
    std::string ret;

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (S[i] == T[j])
            {
                if (i == 0 || j == 0)
                    L[i][j] = 1;
                else
                    L[i][j] = L[i - 1][j - 1] + 1;
 
                if (L[i][j] > z)
                {
                    z = L[i][j];
                    ret = S.substr(i - z + 1, z);
                }
            }
            else
            {
                L[i][j] = 0;
            }
        }
    }
    return ret;
} 
')




mypath <- "/home/owner/Nextcloud/laptop_rplace/eichele/art_marzooplclk/ncbi_submission"
transdir <- paste0(mypath, "/", "transcriptomes")

samps <- list.files(path = transdir, pattern = "[A-Za-z]+_[A-Za-z0-9]+", 
                    all.files = TRUE, full.names = TRUE, include.dirs = TRUE)

#VecScreen can end up going multiple rounds.
#So if there are new errors, just add them at the end of errors1.txt
#And rerun this script.

for(i in 1:length(samps)){
  
  #  
  i <- 3
  cursamp <- samps[[i]]
  cat("Processing errors for", basename(cursamp), "\n")
  
  seqdat <- fastodf(paste0(cursamp, "/", basename(cursamp), "_transcriptome.fsa"))
  errdat <- data.table::fread(paste0(cursamp, "/", "errors1.txt"), sep = "$", header = FALSE)
  
  if(!purrr::is_empty(errdat)){
    
    #Basic munging to get a uid column.
    seqdat %<>% mutate(filename = str_replace(filename, "_transcriptome.fsa", ""))
    seqdat %<>% mutate(uid = paste0(filename, "__", curname))
    
    #Basic munging to get a uid column, and extract
    #the sequence ID + interval from the VecScreen report.
    errdat$samp <- paste0(cursamp, "/", "errors1.txt")
    errdat %<>% mutate(curname = str_extract(V1, "(?<=Sequence-id: ).*?(?=,)"),
                       interval = str_extract(V1, "(?<=Interval: ).*?(?=,)"), 
                       filename = basename(dirname(samp)))
    errdat %<>% select(c(filename, curname, interval))
    errdat %<>% 
      mutate(int_start = as.numeric(str_extract(interval, "^\\d+")), 
             int_stop = as.numeric(str_extract(interval, "\\d+$"))) %>%
      select(-interval)
    errdat %<>% mutate(uid = paste0(filename, "__", curname), .keep = "unused")
    
    
    #Merging the interval data into the sequence data.frame.
    seqdat <- full_join(seqdat, errdat, by = "uid")
    rm(errdat)
    
    #There can be cases where more than one modification has to be
    #performed upon the same sequence. So I will add row_number()
    #to the uid to distinguish these.
    seqdat %<>% mutate(uid = paste0(uid, "__", row_number()))
    
    #Getting the length of the sequence as a column.
    seqdat %<>% mutate(slen = nchar(curseq))
    
    #What is the action to be performed in terms of
    #the contaminant sequence.
    #If no interval is defined, do nothing.
    #If the interval is internal, drop the sequence entirely.
    #If the interval is flanking, trim the sequence.
    seqdat %<>% mutate(act = case_when(
      is.na(int_start) ~ "nothing",
      int_start > 1 & int_stop < slen ~ "drop",
      int_start == 1 | int_stop == slen ~ "trim"
    ))
    
    
    
    #Put the sequences with no changes to be made in a separate df.
    nomoddat <- seqdat %>% 
      filter(act == "nothing") %>% 
      transmute(curseq = paste0(">", curname, "\n", curseq))
    #Prepping this for writing out already, to save space in memory.
    
    #Retain sequences to be modified in seqdat.
    seqdat %<>% filter(act != "nothing")
    
    
    
    
    #If the sequence is internal, we can just extend
    #the trimming window to the nearest end on the parent sequence.
    #So the only action becomes trimming (in the case of a
    #contaminant subsequence).
    seqdat %<>% mutate(dbeg = int_start - 1, 
                       dend = slen - int_stop, 
                       to_trim = ifelse(dbeg > dend, "dend", "dbeg"))
    
    #Pivoting longer as this makes this easier.
    seqdat %<>% pivot_longer(cols = c(int_start, int_stop), names_to = "interval", values_to = "int_pos")
    #Basically look at which end has to be trimmed, and check
    #which interval position this corresponds to (start or stop
    #of the interval), and update the positional data in a new
    #column accordingly.
    seqdat %<>% 
      rowwise() %>%
      mutate(new_int = case_when(
      to_trim == "dbeg" & interval == "int_start" ~ 1,
      to_trim == "dend" & interval == "int_stop" ~ as.numeric(slen),
      TRUE ~ int_pos
    )) %>%
      data.frame(stringsAsFactors = FALSE)
    #At this point, one can check with the interval, int_pos, and new_int
    #columns to see if the trimming window has been extended to the correct
    #flank.
    
    #tst <- seqdat
    #seqdat <- tst
    
    #Removing unused columns.
    seqdat %<>% select(-c(dbeg, dend, to_trim, act, int_pos))
    
    
    #Pivoting wider to get one row per sequence again.
    seqdat %<>% pivot_wider(names_from = interval, values_from = new_int)
    
    #We only need the coordinates of the subsequence to be retained.
    #So converting to that.
    cat("Calculating trimming regions.\n")
    seqdat %<>% 
      rowwise() %>%
      mutate(ret_start = ifelse(int_start == 1, int_stop + 1, 1),
             ret_stop = ifelse(int_stop == slen, int_start - 1, slen), 
             ) %>% 
      data.frame(stringsAsFactors = FALSE)
    
    cat("Trimming.\n")
    seqdat %<>% 
      rowwise() %>%
      mutate(newseq = case_when(
      !is.na(int_start) ~ str_sub(curseq, start = ret_start, end = ret_stop),
      TRUE ~ curseq
    )) %>% 
      data.frame(stringsAsFactors = FALSE)
    
    seqdat %<>% mutate(newlen = nchar(newseq))
    seqdat %<>% mutate(diff = abs(newlen - slen))
    
    
    #seqdat %<>% select(-c(uid, curseq))
    
    cat("Trimmed sequences are:\n")
    trimdat <- seqdat %>% filter(diff > 0) %>% select(-c(newseq, filename, curseq, curname)) %>% print()
    #seqdat %<>% filter(act != "nothing")
    
    #Accumulating trimdat data to a data.frame for
    #later perusal.
    if(!exists("outdf")){
      outdf <- trimdat
    } else{
      outdf <- bind_rows(outdf, trimdat)
    }
    
    
    #If there are multiple edits to the same sequence,
    #the results need to be compared and the
    #longest common substring between the two edits needs to be
    #taken as the representative sequence.
    #If not, the sequence is discarded.
    seqdat %<>% select(c(uid, curname, newseq, slen, int_start, int_stop, ret_start, ret_stop))
    
    #First grouping by sequence ID and collapsing multiple sequences into a single row
    #separated by "__"
    cat("Extracting common substrings for multiple trimming operations, if any exist.\n")
    seqdat %<>%
      group_by(curname) %>%
      mutate(filtseq = paste0(newseq, collapse = "__")) %>%
      ungroup()
    
    #Calling the custom largest_common_substring() function
    #on the two sequences separated by "__".
    #If there aren't multiple trimmed versions of the strings
    #to be compared, nothing is done.
    seqdat %<>% 
      rowwise() %>%
      mutate(filtseq = ifelse(str_detect(filtseq, "__"), 
                                       largest_common_substring(unlist(str_split(filtseq, "__"))[[1]], 
                                                                unlist(str_split(filtseq, "__"))[[1]]), 
                              filtseq)) %>%
      data.frame() %>%
      select(-c(newseq, int_start, int_stop, ret_start, ret_stop, slen))
    #Also compacting the data.frame here.
    
    #Retaining just one copy of the longest (final) common trimmed
    #sequence (since each row would have copy).
    #If there are any sequences that are SHORTER than 200 bp
    #after trimming, they will be discarded (NCBI limit).
    seqdat %<>% 
      group_by(curname) %>%
      distinct(filtseq, .keep_all = TRUE) %>%
      ungroup() %>%
      filter(nchar(filtseq) >= 200)
    
    
    #Merging seqdat and nomoddat together.
    seqdat %<>% transmute(curseq = paste0(">", curname, "\n", filtseq))
    seqdat %<>% bind_rows(nomoddat)
    rm(nomoddat)
    
    #Writing out the data.
    sname <- paste0(cursamp, "/", basename(cursamp), "_transcriptome_cleaned.fsa")
    #seqdat %<>% transmute(newseq = paste0(">", curname, "\n", newseq))
    
    cat("Writing out", sname, "\n")
    write.table(seqdat, file = sname, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    cat("Done.\n")
  } else{
    #If there were no VecScreen errors.
    #Write out the input as the cleaned fasta file.
    seqdat %<>% transmute(curseq = paste0(">", curname, "\n", curseq))
    
    cat("No errors from VecScreen. Skipping.\n")
    sname <- paste0(cursamp, "/", basename(cursamp), "_transcriptome_cleaned.fsa")
    #seqdat %<>% transmute(newseq = paste0(">", curname, "\n", newseq))
    
    cat("Writing out", sname, "\n")
    write.table(seqdat, file = sname, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    cat("Done.\n")
    
  }
  
}

cat("All operations done.\n")

#tbl2asn execution commandline call:
#Call from the transcriptomes directory.
#Make sure the the tbl2asn_conda environment is activated.
#for FLE in ./*; do SAMP=$(basename $FLE); echo $SAMP; SAMPINS=${SAMP/_/ }; SAMPINS=$(echo ${SAMPINS} | sed -r "s/spHR/sp\. HR/"); cd ${FLE}; echo ${PWD}; tbl2asn -p. -t template.sbt -w assembly.cmt -M t -i ${SAMP}_transcriptome_cleaned.fsa -j "[organism=${SAMPINS}][moltype=mRNA][tech=TSA]" -V t -W; cd ..; echo ${PWD}; done