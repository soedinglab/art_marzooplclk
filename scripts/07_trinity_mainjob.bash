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
TRINJOBDIR="${MAINJOBDIR}/trinity/jobsubs";
mkdir -p ${TRINJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
TRINOUTDIR="${MAINOUTDIR}/trinity";
mkdir -p ${TRINOUTDIR};

NVMEPATH="/nvme/n00/user0/newrun";

#READSPATH="${MYPATH}/reads";
READSPATH="${MAINOUTDIR}/sortmerna";


for SAMP in $(ls ${READSPATH}/*.fastq | sed -r 's/_R[12].*\.fastq$//' | uniq); 
do

    R1=${SAMP}_R1_smr.fastq;
    R2=${SAMP}_R2_smr.fastq;
    BR1=$(basename ${R1});
    BR2=$(basename ${R2});
    TOUTDIR=$(basename ${TRINOUTDIR});

    #R1_FP=$(echo $(basename ${R1}) | sed 's/_kra2/_fp/');
    #R2_FP=$(echo $(basename ${R1}) | sed 's/_kra2/_fp/');

    SAMPBNAME=$(basename ${SAMP});
    
    JOB_FILE="${TRINJOBDIR}/${SAMPBNAME}_trinity_jobsub.bash";

    #TRINSING="/disk0/user0/my_programs/singularity_images/trinityrnaseq.v2.11.0.simg";
    TRINSING="/disk0/user0/my_programs/singularity_images/trinityrnaseq.v2.12.0.simg";

    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J_${SAMPBNAME}
#SBATCH -e err.minijob.%J_${SAMPBNAME}
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH -C "rome"
#SBATCH --exclude=hh001
#SBATCH --hint=nomultithread


#mkdir -p ${TRINOUTDIR}/${SAMPBNAME}_trinity;
#cd ${TRINOUTDIR}/${SAMPBNAME}_trinity;

#mkdir -p "${NVMEPATH}/reads";
#cp ${R1} "${NVMEPATH}/reads";
#cp ${R2} "${NVMEPATH}/reads";

#mkdir -p ${NVMEPATH}/outputs/${TOUTDIR}/${SAMPBNAME}_trinity;
#cd ${NVMEPATH}/outputs/${TOUTDIR}/${SAMPBNAME}_trinity;


#
source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;

#module load singularity/3.7.0;
#singularity exec -B ${NVMEPATH} -e ${TRINSING} Trinity --seqType fq --left ${NVMEPATH}/reads/${BR1} --right ${NVMEPATH}/reads/${BR2} --output ${NVMEPATH}/outputs/${TOUTDIR}/${SAMPBNAME}_trinity --full_cleanup --verbose --max_memory 240G --bflyCalculateCPU --CPU 32;
#cp ${NVMEPATH}/outputs/${TOUTDIR}/${SAMPBNAME}_trinity.Trinity.fasta ${TRINOUTDIR};

#Need a minimum length of 200 for TSA submissions, but Trinity isn't playing nice in that regard, so filtering for that.
conda activate seqkit_conda;
seqkit seq -w 0 -m 200 -j 32 -o ${TRINOUTDIR}/${SAMPBNAME}_trinity.fasta ${TRINOUTDIR}/${SAMPBNAME}_trinity.Trinity.fasta;
conda deactivate;

    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done

