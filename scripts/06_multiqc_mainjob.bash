#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 14-00:00:00
#SBATCH -C "haswell"


#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
MQCJOBDIR="${MAINJOBDIR}/multiqc/jobsubs";
mkdir -p ${MQCJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
MQCOUTDIR="${MAINOUTDIR}/multiqc";
mkdir -p ${MQCOUTDIR};


for FQCDIR in ${MAINOUTDIR}/fastqc_*; 
do

    SAMPBNAME=$(basename ${FQCDIR});

    source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
    
    conda activate multiqc_conda;

    multiqc --no-data-dir --flat --fullnames --outdir ${MQCOUTDIR} --title ${SAMPBNAME} ${FQCDIR};
    
    conda deactivate;

done
