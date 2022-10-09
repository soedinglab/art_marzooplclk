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
MPCLUJOBDIR="${MAINJOBDIR}/mmseqs2_protclu_td/jobsubs";
mkdir -p ${MPCLUJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
MPCLUOUTDIR="${MAINOUTDIR}/mmseqs2_protclu_td";
mkdir -p ${MPCLUOUTDIR};


#READSPATH="${MYPATH}/reads";
FILESPATH="${MAINOUTDIR}/transdecoder_trinity";


for SAMPFA in ${FILESPATH}/*_transdecoder_fixed_noast.fasta; 
do

    SAMPBNAME=$(basename ${SAMPFA} _transdecoder_fixed_noast.fasta);

    JOB_FILE="${MPCLUJOBDIR}/${SAMPBNAME}_mmseqs2_protclu_trinity_jobsub.bash";



    #Two stage clustering
    CLUSAMP1=${SAMPBNAME}_pclu;
    #CLUSAMP2=${SAMPBNAME}_pclu2;

    #tmp directories
    CLUTMP1="${MPCLUOUTDIR}/${SAMPBNAME}_pclu_tmp";
    #CLUTMP2="${MPCLUOUTDIR}/${SAMPBNAME}_pclu2_tmp";
    mkdir -p ${CLUTMP1};
    #mkdir -p ${CLUTMP2};
    
    #Rep seqs from each round
    CLUSAMP1REPSEQ=${MPCLUOUTDIR}/${CLUSAMP1}_rep_seq.fasta;
    #CLUSAMP2REPSEQ=${MPCLUOUTDIR}/${CLUSAMP2}_rep_seq.fasta;


    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=124G
#SBATCH --hint=nomultithread
#SBATCH -C "haswell"

source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
conda activate mmseqs2_conda;

cd ${MPCLUOUTDIR};

#Clustering round 1
mmseqs easy-linclust ${SAMPFA} ${CLUSAMP1} ${CLUTMP1} --kmer-per-seq 100 --kmer-per-seq-scale 0.8 --comp-bias-corr 0 --mask 0 -e 1000000000 --cov-mode 1 --cluster-mode 2 -c 0.98 --min-seq-id 0.98;

#Clustering round 2
#mmseqs easy-linclust ${CLUSAMP1REPSEQ} ${CLUSAMP2} ${CLUTMP2} --kmer-per-seq 100 --kmer-per-seq-scale 0.8 --comp-bias-corr 0 --mask 0 -e 1000000000 --cov-mode 1 --cluster-mode 2 -c 0.90 --min-seq-id 0.98

conda deactivate;


    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done
