#!/bin/bash

#Paths and stuff
MYPATH="/disk0/user0/newrun";

TNAME="rscript_ncbicomp";

MAINJOBDIR="${MYPATH}/jobfiles";
TASKJOBDIR="${MAINJOBDIR}/${TNAME}/jobsubs";
mkdir -p ${TASKJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
TASKOUTDIR="${MAINOUTDIR}/${TNAME}";
mkdir -p ${TASKOUTDIR};


#Species list.
SPLIST="${MYPATH}/species_taxids.csv";
#Circadian clock proteins queries fasta file.
QFA="${MYPATH}/cc_queries.fasta";

#Sourcing conda.
source /path/to/miniconda3/etc/profile.d/conda.sh;


#The -u4 here and the 4< to feed into th
#while loop are crucial.
#See https://stackoverflow.com/a/66179139 for
#an explanation.
while IFS="," read -u4 ORG TXID CMT;
do
    
    SNAME=$(echo ${ORG/ /_});
    SNAME=$(echo ${SNAME/./});

    #SOUTDIR="${TASKOUTDIR}/${SNAME}";
    #mkdir -p ${SOUTDIR};

    echo ${SNAME};
    
    conda activate entrez_conda;
    #Search with entrez using the TXIDs.
    QUERY='((txid'${TXID}'[Organism:exp]) AND ( "tsa-master"[Properties] OR "wgs-master"[Properties] ))';
    echo $QUERY;
    RESULT=$(esearch -db nuccore -query "${QUERY}" | efetch -format xml | tee ${TASKOUTDIR}/${SNAME}_${TXID}.esearch.xml);
    #echo ${RESULT};
    # save the xml result for further reference
    IDS=$(echo $RESULT | xtract -pattern Seq-entry -element Textseq-id_name);
    conda deactivate;

    if [ -n "${IDS}" ]
    then
        echo "Something available, running fastq-dump!!";
        conda activate entrez_conda;
        for ID in ${IDS}; 
        do
            FA="${TASKOUTDIR}/${SNAME}_ncbi_${ID}.fasta";
            if [ ! -f ${FA} ]
            then
                cd ${TASKOUTDIR};
                # use SRA toolkit fastq-dump with the option to make a fasta file and the standard header.
                echo "Downloading data for ${FA}!!";
                fastq-dump --fasta 0 -F ${ID};
                mv ${TASKOUTDIR}/${ID}.fasta ${FA};
                echo "Done downloading data for ${FA}!!";
            else
                echo "File already exists, skipping download ${FA}!!";
            fi
        done
        conda deactivate;
    else
        echo "Nothing available";
    fi
    echo "Inner done!!";
done 4< ${SPLIST};


#MMseqs2 search
conda activate mmseqs2_conda
for TFA in ${TASKOUTDIR}/*.fasta;
do
    OUTTAB=${TFA/.fasta/.tab};
    if [ ! -f ${OUTTAB} ]
    then
        echo "Searching against ${TFA}!!";
        mmseqs easy-search ${QFA} ${TFA} ${OUTTAB} ${OUTTAB}_tmp -s 7.5 -e 0.00001 --format-output qheader,theader,pident,evalue,bits,qlen,tlen,alnlen,qcov,tcov;
        echo "Done searching against ${TFA}!!";
    else
        echo "Search results exist already: ${OUTTAB}!!"
    fi
done
conda deactivate;

#Putting the downloaded 
ls ${TASKOUTDIR}/*_ncbi_*.fasta | tee ${TASKOUTDIR}/"ncbi_downloaded_genomes_proteomes.txt";
echo "Outer done!!";
