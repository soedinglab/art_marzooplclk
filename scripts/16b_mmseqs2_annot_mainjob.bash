#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 14-00:00:00
#SBATCH -C "haswell"



#Reference annotation database.
#Just the UniProt SwissProt database.
TDSP="/disk0/user0/my_programs/transdecoder_dbs/transdecoder_dbs_diamond/uniprot_sprot.dmnd";
TDSPFA="/disk0/user0/my_programs/transdecoder_dbs/transdecoder_dbs_diamond/uniprot_sprot.fasta";

#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
TDJOBDIR="${MAINJOBDIR}/mmseqs2_annots_rfilt/jobsubs";
mkdir -p ${TDJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
MAINTDOUTDIR="${MAINOUTDIR}/mmseqs2_annots_rfilt";
mkdir -p ${MAINTDOUTDIR};

#Local outdir. This is what $SAMPTDOUTDIR is under.
NVMEPATH="/local_hdd/user0/newrun";
TDOUTDIR="${NVMEPATH}/outputs/mmseqs2_annots_rfilt";
mkdir -p ${TDOUTDIR};


FILESPATH="${MAINOUTDIR}/rscript_lenfilt";

for SAMPFA in ${FILESPATH}/*_td_rfilt.fasta; 
do

    SAMPBNAME=$(basename ${SAMPFA} _td_rfilt.fasta);

    JOB_FILE="${TDJOBDIR}/${SAMPBNAME}_mmseqs2_annots_rfilt_jobsub.bash";

    SAMPTDOUTDIR=${TDOUTDIR}/${SAMPBNAME};
    mkdir -p ${SAMPTDOUTDIR};


    echo "#!/bin/bash -ex
#SBATCH -J minijob_td
#SBATCH -p hh
#SBATCH -o out.minijob.%J_${SAMPBNAME}
#SBATCH -e err.minijob.%J_${SAMPBNAME}
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=124G
#SBATCH --hint=nomultithread
#SBATCH -C "haswell"

mkdir -p ${SAMPTDOUTDIR};
cd ${SAMPTDOUTDIR};


#Local copies of everything I need.
cp ${SAMPFA} ${SAMPTDOUTDIR}/$(basename ${SAMPFA});
cp ${TDSP} ${SAMPTDOUTDIR}/$(basename ${TDSP});
cp ${TDSPFA} ${SAMPTDOUTDIR}/$(basename ${TDSPFA});


source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;


#Note supplying outfmt6 with custom headers
conda activate mmseqs2_conda;
mmseqs easy-search ${SAMPTDOUTDIR}/$(basename ${SAMPFA}) ${SAMPTDOUTDIR}/$(basename ${TDSPFA}) ${SAMPTDOUTDIR}/${SAMPBNAME}_mmseqs2_outfmt6.tab ${SAMPTDOUTDIR}/${SAMPBNAME}_tmp -s 7.5 -e 0.00001 --format-output qheader,theader,pident,evalue,bits,qlen,tlen,alnlen,qcov,tcov;
conda deactivate;

conda activate diamond_conda;
diamond blastp --ultra-sensitive --query ${SAMPTDOUTDIR}/$(basename ${SAMPFA}) --db ${SAMPTDOUTDIR}/$(basename ${TDSP}) --max-target-seqs 1 --outfmt 6 qtitle stitle pident evalue bitscore qlen slen length --header --evalue 1e-5 --threads 16 --out ${SAMPTDOUTDIR}/${SAMPBNAME}_dmnd_outfmt6.tab;
conda deactivate;


#Copying over the annotations table to /disk0.
cp ${SAMPTDOUTDIR}/${SAMPBNAME}_mmseqs2_outfmt6.tab ${MAINTDOUTDIR};
cp ${SAMPTDOUTDIR}/${SAMPBNAME}_dmnd_outfmt6.tab ${MAINTDOUTDIR};

    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done