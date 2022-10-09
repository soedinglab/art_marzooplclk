#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00
#SBATCH --hint=nomultithread


#module load JAVA/OPENJDK/13.0.1
#fastqc=$(realpath $HOME)"/my_programs/FastQC/fastqc";

#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
KRAJOBDIR="${MAINJOBDIR}/kraken2/jobsubs";
mkdir -p ${KRAJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
KRAOUTDIR="${MAINOUTDIR}/kraken2";
mkdir -p ${KRAOUTDIR};

#Kraken2 database
KRADB="/disk0/user0/my_programs/kraken2_dbs/pluspfp/";
KRATOOLS="/usr/users/user0/my_programs/KrakenTools/";
EXTKRAREADS="${KRATOOLS}/extract_kraken_reads.py"

#READSPATH="${MYPATH}/reads";
READSPATH="${MAINOUTDIR}/fastp";


for SAMP in $(ls ${READSPATH}/*.fastq | sed -r 's/_R[12].*\.fastq$//' | uniq); 
do

    R1=${SAMP}_R1_fp.fastq;
    R2=${SAMP}_R2_fp.fastq;
    
    JOB_FILE="${KRAJOBDIR}/$(basename ${SAMP})_fastqc_kra_jobsub.bash";

    SAMPBNAME=$(basename ${SAMP});

    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=900G
#SBATCH --exclude=hh001

#source /usr/users/user0/miniconda3/etc/profile.d/conda.sh
source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
conda activate kraken2_conda;

mkdir -p ${KRAOUTDIR}/${SAMPBNAME};
cd ${KRAOUTDIR}/${SAMPBNAME};

#Classified outs = contaminants
COUT="${SAMPBNAME}_kra2_contaminants_R";
#Unclassified outs = clean reads
UCOUT="${SAMPBNAME}_kra2_R";

#Kraken2 needs to hold the entire DB in RAM, so I will give it a ton of memory and threads
#If I use --memory-mapping then the DB would not be held in RAM but this is much slower
kraken2 --minimum-hit-groups 4 --confidence 0.10 --report-zero-counts --use-names --paired --threads 128 --db ${KRADB} --report ${SAMPBNAME}_kraken2_report --classified-out ${SAMPBNAME}_R#_kra2_cont.fastq --unclassified-out ${SAMPBNAME}_R#_kra2.fastq ${R1} ${R2}

conda deactivate;

#
mv ${KRAOUTDIR}/${SAMPBNAME}/* ${KRAOUTDIR}

    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done

#for FLE in ./*.fastq; do echo mv ${FLE} $(dirname ${FLE})/$(echo $(basename ${FLE}) | sed -r 's/_R_/_R/'); done
#Manual bash one-liner to rename files properly at the end.
