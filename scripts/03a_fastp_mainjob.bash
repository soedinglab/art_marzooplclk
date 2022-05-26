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
FPJOBDIR="${MAINJOBDIR}/fastp/jobsubs";
mkdir -p ${FPJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
FPOUTDIR="${MAINOUTDIR}/fastp";
mkdir -p ${FPOUTDIR};


#READSPATH="${MYPATH}/reads";
#
READSPATH="${MAINOUTDIR}/rcorrector";


for SAMP in $(ls ${READSPATH}/*.fastq | sed -r 's/_R[12].*\.fastq$//' | uniq); 
do

    R1=${SAMP}_R1_rcor.fastq;
    R2=${SAMP}_R2_rcor.fastq;

    R1_FP=$(basename ${R1} _rcor.fastq)_fp.fastq;
    R2_FP=$(basename ${R2} _rcor.fastq)_fp.fastq;
    
    JOB_FILE="${FPJOBDIR}/$(basename ${SAMP})_fastp_jobsub.bash";

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

source /usr/users/user0/miniconda3/etc/profile.d/conda.sh
conda activate fastp_conda

#mkdir -p ${FPOUTDIR}/${SAMPBNAME};
#cd ${FPOUTDIR}/${SAMPBNAME};

cd ${FPOUTDIR};

fastp -i ${R1} -I ${R2} -o ${R1_FP} -O ${R2_FP} --json "${SAMPBNAME}_fastp.json" --html "${SAMPBNAME}_fastp.html" --thread=16 --detect_adapter_for_pe --correction --overrepresentation_analysis --low_complexity_filter --trim_poly_x --trim_poly_g --average_qual=20 --qualified_quality_phred=20 --overrepresentation_sampling=1 --n_base_limit 0 --length_required 35;

conda deactivate


    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done
