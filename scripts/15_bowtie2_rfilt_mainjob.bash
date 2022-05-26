#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00
#SBATCH --hint=nomultithread


#module load JAVA/OPENJDK/13.0.1


#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
BT2JOBDIR="${MAINJOBDIR}/bowtie2_rfilt/jobsubs";
mkdir -p ${BT2JOBDIR};


MAINOUTDIR="${MYPATH}/outputs";
MAINBT2OUTDIR="${MAINOUTDIR}/bowtie2_rfilt";
mkdir -p ${MAINBT2OUTDIR};

NVMEPATH="/local_hdd/user0/newrun";
BT2OUTDIR="${NVMEPATH}/outputs/bowtie2_rfilt";
#mkdir -p ${BT2OUTDIR};


#READSPATH="${MYPATH}/reads";
READSPATH="${MAINOUTDIR}/sortmerna";

#Assemblies path -- these are the transcriptomes filtered
#on the basis of rfilt proteomes.
ASSEMPATH="${MAINOUTDIR}/rscript_annotscmp/transcriptomes";


for SAMPFA in ${ASSEMPATH}/*_transcriptome.fasta; 
do

    SAMPBNAME=$(basename ${SAMPFA} _transcriptome.fasta);

    #File name shenanigans since the rfilt transcriptomes do not have the grp suffix in them.
    if [ "${SAMPBNAME}" == "Acartia_tonsa" ] || [ "${SAMPBNAME}" == "Acartia_clausii" ] || [ "${SAMPBNAME}" == "Evadne_nordmanni" ] || [ "${SAMPBNAME}" == "Oikopleura_dioica" ] || [ "${SAMPBNAME}" == "Temora_longicornis" ];
    then
        SAMPBNAME="${SAMPBNAME}grp";
    fi;

    R1=${READSPATH}/${SAMPBNAME}_R1_smr.fastq;
    R2=${READSPATH}/${SAMPBNAME}_R2_smr.fastq;



    JOB_FILE="${BT2JOBDIR}/${SAMPBNAME}_bowtie2_jobsub.bash";

    SAMPBT2OUTDIR=${BT2OUTDIR}/${SAMPBNAME};
    BR1=${SAMPBT2OUTDIR}/$(basename ${R1});
    BR2=${SAMPBT2OUTDIR}/$(basename ${R2});
    LOCSAMPFA=${SAMPBT2OUTDIR}/$(basename ${SAMPFA});



    echo "#!/bin/bash -ex
#SBATCH -J minijob_bt2
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=124G
#SBATCH --hint=nomultithread
#SBATCH -C "haswell"


source /usr/users/user0/miniconda3/etc/profile.d/conda.sh;
conda activate bowtie2_conda;

mkdir -p ${SAMPBT2OUTDIR};
cd ${SAMPBT2OUTDIR};

cp ${R1} ${BR1};
cp ${R2} ${BR2};
cp ${SAMPFA} ${LOCSAMPFA};

#Building bowtie2 index
bowtie2-build ${LOCSAMPFA} ${SAMPBT2OUTDIR}/${SAMPBNAME};

#Aligning the reads
#Basically following this https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly
bowtie2 --no-unal --threads 16 -q -k 20 -x ${SAMPBT2OUTDIR}/${SAMPBNAME} -1 ${BR1} -2 ${BR2} 2>${SAMPBT2OUTDIR}/${SAMPBNAME}_align_stats.txt | samtools view -@10 -Sb -o ${SAMPBT2OUTDIR}/${SAMPBNAME}_bowtie2.bam;


conda deactivate;


mkdir -p ${MAINBT2OUTDIR}/${SAMPBNAME};
cp ${SAMPBT2OUTDIR}/${SAMPBNAME}_bowtie2.bam ${MAINBT2OUTDIR}/${SAMPBNAME};
cp ${SAMPBT2OUTDIR}/*.bt2 ${MAINBT2OUTDIR}/${SAMPBNAME};
cp ${SAMPBT2OUTDIR}/*.txt ${MAINBT2OUTDIR}/${SAMPBNAME};



    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done
