#!/bin/bash -ex
#SBATCH -J mainjob
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00


#Paths and stuff
MYPATH="/disk0/user0/newrun";


#Job files paths
MAINJOBDIR="${MYPATH}/jobfiles";
SMRJOBDIR="${MAINJOBDIR}/sortmerna/jobsubs";
mkdir -p ${SMRJOBDIR};

#Sortmerna DBs
#SR_DB="/home/somedir/user0/my_programs/sortmerna/data/rRNA_databases";
SMR_DBP="/disk0/user0/jobfiles/sortmerna/rnadbs";


#DB1="--ref ${SMR_DBP}/rfam-5.8s-database-id98.fasta";
#DB2="--ref ${SMR_DBP}/rfam-5s-database-id98.fasta";
#DB3="--ref ${SMR_DBP}/silva-arc-16s-id95.fasta";
#DB4="--ref ${SMR_DBP}/silva-arc-23s-id98.fasta";
#DB5="--ref ${SMR_DBP}/silva-bac-16s-id90.fasta";
#DB6="--ref ${SMR_DBP}/silva-bac-23s-id98.fasta";
#DB7="--ref ${SMR_DBP}/silva-euk-18s-id95.fasta";
#DB8="--ref ${SMR_DBP}/silva-euk-28s-id98.fasta";

#DBS="${DB1} ${DB2} ${DB3} ${DB4} ${DB5} ${DB6} ${DB7} ${DB8}";

DBS="--ref /disk0/user0/my_programs/sortmerna/data/smr_v4.3_sensitive_db.fasta";


#Output paths
MAINOUTDIR="${MYPATH}/outputs";
SMROUTDIR="${MAINOUTDIR}/sortmerna";
mkdir -p ${SMROUTDIR};


#Reads for input
READSPATH="${MAINOUTDIR}/kraken2";

for SAMP in $(ls ${READSPATH}/*.fastq | sed -r 's/_R[12].*.fastq$//' | uniq); 
do

    R1=${SAMP}_R1_kra2.fastq;
    R2=${SAMP}_R2_kra2.fastq;

    SAMPBNAME=$(basename ${SAMP});
    
    JOB_FILE=${SMRJOBDIR}/${SAMPBNAME}_sortmerna_jobsub.bash;
    
    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=900G
#SBATCH --hint=nomultithread
#SBATCH --exclude=hh001
#SBATCH -C "rome"

source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
conda activate sortmerna_conda

mkdir -p ${SMROUTDIR}/${SAMPBNAME};
cd ${SMROUTDIR}/${SAMPBNAME};

#Chose --paired_in based on text after table 3
#here https://github.com/biocore/sortmerna/wiki/User-manual-v4.0#filtering-paired-end-reads
#
#/home/somedir/user0/my_programs/sortmerna-4.2.0-Linux/bin/sortmerna
#
sortmerna --reads ${R1} --reads ${R2} --workdir ${SMROUTDIR}/${SAMPBNAME} ${DBS} --threads 128 --fastx True --paired_in True --out2 True --other;

#File renaming + moving
mv ${SMROUTDIR}/${SAMPBNAME}/out/other_fwd.fq ${SMROUTDIR}/$(basename ${R1} _kra2.fastq)_smr.fastq;
mv ${SMROUTDIR}/${SAMPBNAME}/out/other_rev.fq ${SMROUTDIR}/$(basename ${R2} _kra2.fastq)_smr.fastq;
mv ${SMROUTDIR}/${SAMPBNAME}/out/aligned_fwd.fq ${SMROUTDIR}/$(basename ${R1} _kra2.fastq)_smr_rrna.fastq;
mv ${SMROUTDIR}/${SAMPBNAME}/out/aligned_rev.fq ${SMROUTDIR}/$(basename ${R2} _kra2.fastq)_smr_rrna.fastq;

mv ${SMROUTDIR}/${SAMPBNAME}/out/aligned.log ${SMROUTDIR}/${SAMPBNAME}_aligned.log;
#rm -drf ${SMROUTDIR}/${SAMPBNAME};

conda deactivate

    " > ${JOB_FILE};
  sbatch ${JOB_FILE};

done
