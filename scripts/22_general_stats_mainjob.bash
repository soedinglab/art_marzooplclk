#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH --cpus-per-task=256
#SBATCH --mem=900G
#SBATCH -t 14-00:00:00
#SBATCH -C "rome"

#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
GSJOBDIR="${MAINJOBDIR}/general_stats/jobsubs";
mkdir -p ${GSJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
GSOUTDIR="${MAINOUTDIR}/general_stats";
mkdir -p ${GSOUTDIR};


source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
conda activate seqkit_conda;



#BowTie2 statistics.
echo "Getting Bowtie2 statistics!!\n";
find ${MAINOUTDIR}/bowtie2_* -maxdepth 2 -type f -name "*_stats.txt" -printf "%p " -exec grep "overall alignment rate$" {} \; > ${GSOUTDIR}/bowtie2_align_stats.tsv;
#| sed -r 's/.*(bowtie2_[a-z]+)\/([A-Za-z]+_[a-z]+).*\s([0-9]+\.[0-9]+).*$/\2\t\1\t\3/g' 



#BUSCO statistics.
echo "Getting BUSCO statistics!!\n";
find ${MAINOUTDIR}/busco_* -maxdepth 2 -type f -name "short_summary*txt" -printf "%p " -exec grep "C:" {} \; > ${GSOUTDIR}/busco_all_stats.tsv;

#Read counts.
echo "Getting read counts!!\n";

#Raw reads.
echo "Raw reads!!\n";
seqkit stats -b -T -j 256 --out-file ${GSOUTDIR}/seqkit_raw_stats.tsv ${MYPATH}/reads/*.fastq;

#Rcorrector reads.
echo "Rcorrector reads!!\n";
seqkit stats -b -T -j 256 --out-file ${GSOUTDIR}/seqkit_rcor_stats.tsv ${MAINOUTDIR}/rcorrector/*_rcor.fastq;

#Fastp reads.
echo "Fastp reads!!\n";
seqkit stats -b -T -j 256 --out-file ${GSOUTDIR}/seqkit_fp_stats.tsv ${MAINOUTDIR}/fastp/*_fp.fastq;

#Kraken2 reads.
echo "Kraken2 reads!!\n";
seqkit stats -b -T -j 256 --out-file ${GSOUTDIR}/seqkit_kra2_stats.tsv ${MAINOUTDIR}/kraken2/*_kra2.fastq;

#SortMeRNA reads.
echo "SortMeRNA reads!!\n";
seqkit stats -b -T -j 256 --out-file ${GSOUTDIR}/seqkit_smr_stats.tsv ${MAINOUTDIR}/sortmerna/*_smr.fastq;

echo "Done!!\n";


conda deactivate;

