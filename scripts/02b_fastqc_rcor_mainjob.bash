#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00
#SBATCH --hint=nomultithread


#module load JAVA/OPENJDK/13.0.1
fastqc=$(realpath $HOME)"/my_programs/FastQC/fastqc";

#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
FQCJOBDIR="${MAINJOBDIR}/fastqc_rcor/jobsubs";
mkdir -p ${FQCJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
FQCOUTDIR="${MAINOUTDIR}/fastqc_rcor";
mkdir -p ${FQCOUTDIR};

#READSPATH="${MYPATH}/reads";
READSPATH="${MAINOUTDIR}/rcorrector";


for SAMP in $(ls ${READSPATH}/*.fastq | sed -r 's/_R[12].*\.fastq$//' | uniq); 
do

    R1=${SAMP}_R1_rcor.fastq;
    R2=${SAMP}_R2_rcor.fastq;
    
    JOB_FILE="${FQCJOBDIR}/$(basename ${SAMP})_fastqc_rcor_jobsub.bash";

    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=7G
#SBATCH --hint=nomultithread
#SBATCH -C "rome"
#SBATCH --exclude=hh007

#source /usr/users/user0/miniconda3/etc/profile.d/conda.sh
#conda activate mmseqs_conda
source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;

conda activate fastqc_conda
fastqc "${R1}" "${R2}" --outdir "${FQCOUTDIR}" --threads 16;
conda deactivate
    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done
