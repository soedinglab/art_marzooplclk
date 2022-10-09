#!/bin/bash

#Fine name fixes: basically get only TRINITY_DNwhatever_*.pX_orftype_*_len_* as the header.
#Do this for both the CDS file and the protein file.
#For protein file, also removing the tailing asterisks if they exist.

#CDS file.
for FLE in $1/*_transdecoder_cds.fasta;
do

echo "Fixing file name for CDS sequence of $(basename ${FLE})!!"; 

sed -r 's/(TRINITY_[A-Za-z0-9\_\.]+).*(ORF\stype\:[A-Za-z0-9\_]+).*/\1_\2/g; s/\s//g; s/\:/_/g' ${FLE} > $(dirname ${FLE})/$(basename ${FLE} .fasta)_fixed.fasta;

echo "Done!! Output in $(basename ${FLE} .fasta)_fixed.fasta!!";

done



#Protein file.
for FLE in $1/*_transdecoder.fasta;
do

echo "Fixing file name for protein sequences of $(basename ${FLE})!!"; 

sed -r 's/(TRINITY_[A-Za-z0-9\_\.]+).*(ORF\stype\:[A-Za-z0-9\_]+).*/\1_\2/g; s/\s//g; s/\:/_/g' ${FLE} > $(dirname ${FLE})/$(basename ${FLE} .fasta)_fixed.fasta;

echo "Replacing all asterisk stop codon characters with nothing in file $(basename ${FLE})!!"; 

sed 's/\*$//g' $(dirname ${FLE})/$(basename ${FLE} .fasta)_fixed.fasta > $(dirname ${FLE})/$(basename ${FLE} .fasta)_fixed_noast.fasta;

echo "Done!! Output in $(basename ${FLE} .fasta)_fixed.fasta!!";

done







