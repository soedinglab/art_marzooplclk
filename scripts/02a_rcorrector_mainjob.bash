#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00
#SBATCH --hint=nomultithread


#module load JAVA/OPENJDK/13.0.1

#Path to Rcorrector post-processing Python script
RCORPOST="/disk0/user0/my_programs/TranscriptomeAssemblyTools/rcorrector_filter_uncor_myver.py";

#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
FPJOBDIR="${MAINJOBDIR}/rcorrector/jobsubs";
mkdir -p ${FPJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
FPOUTDIR="${MAINOUTDIR}/rcorrector";
mkdir -p ${FPOUTDIR};


#
READSPATH="${MYPATH}/reads";
#READSPATH="${MAINOUTDIR}/kraken2";


for SAMP in $(ls ${READSPATH}/*.fastq | sed -r 's/_R[12].*\.fastq$//' | uniq); 
do

    R1=${SAMP}_R1.fastq;
    R2=${SAMP}_R2.fastq;

    R1_FP=$(basename ${R1} .fastq)_fp.fastq;
    R2_FP=$(basename ${R2} .fastq)_fp.fastq;
    
    JOB_FILE="${FPJOBDIR}/$(basename ${SAMP})_rcorrector_jobsub.bash";

    SAMPBNAME=$(basename ${SAMP});

    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=7G
#SBATCH --hint=nomultithread
#SBATCH -C "haswell"

source /usr/users/user0/miniconda3/etc/profile.d/conda.sh;
conda activate rcorrector_conda;

mkdir -p ${FPOUTDIR}/${SAMPBNAME};
cd ${FPOUTDIR}/${SAMPBNAME};

#cd ${FPOUTDIR};
#run_rcorrector.pl -1 ${R1} -2 ${R2} -od ${FPOUTDIR}/${SAMPBNAME} -t 16;

python3 ${RCORPOST} -l ${FPOUTDIR}/${SAMPBNAME}/${SAMPBNAME}_R1.cor.fq -r ${FPOUTDIR}/${SAMPBNAME}/${SAMPBNAME}_R2.cor.fq -o ${FPOUTDIR}

conda deactivate


    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done
