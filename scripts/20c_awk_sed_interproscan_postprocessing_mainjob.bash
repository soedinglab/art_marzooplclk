#!/bin/bash

for FLE in $1/*_interpro_res.tsv;
do
    echo "Fixing TSV file ${FLE} to change blanks to NAs!!";

    awk -F"\t" -v OFS="\t" '{
        for (i=1;i<=15;i++) {
            if ($i == "") $i="NA"
        }
        print $0
    }' ${FLE} > $1/$(basename ${FLE} .tsv)_fixed.tsv;

    echo "Done!!"

done


for FLE in $1/*_interpro_res.gff3;
do

    echo "Fixing GFF3 file ${FLE} by creating FASTA sequence-less version!!";

    sed -n '/^##FASTA.*/q;p' ${FLE} > $1/$(basename ${FLE} .gff3)_fixed.gff3;

    echo "Done!!";


done


