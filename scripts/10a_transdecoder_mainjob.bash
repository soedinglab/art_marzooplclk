#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 14-00:00:00
#SBATCH -C "haswell"


#module load JAVA/OPENJDK/13.0.1


#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
TDJOBDIR="${MAINJOBDIR}/transdecoder_trinity/jobsubs";
mkdir -p ${TDJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
MAINTDOUTDIR="${MAINOUTDIR}/transdecoder_trinity";
mkdir -p ${MAINTDOUTDIR};

#Local outdir. This is what $SAMPTDOUTDIR is under.
NVMEPATH="/local_hdd/user0/newrun";
TDOUTDIR="${NVMEPATH}/outputs/transdecoder_trinity";
mkdir -p ${TDOUTDIR};


#READSPATH="${MYPATH}/reads";
FILESPATH="${MAINOUTDIR}/trinity";
#FILESPATH="${MAINOUTDIR}/infernal_trinity";

#for SAMPFA in ${FILESPATH}/*_trinity.Trinity.fasta; 
#for SAMPFA in ${FILESPATH}/*_inf_mrna.fasta;
#This is the version with the _trinity.Trinity.fasta file with sequences shorter than 200 nucleotides excluded.
#This is necessary for TSA submissions.
for SAMPFA in ${FILESPATH}/*_trinity.fasta;
do

    #Transdecoder paths
    TDPATH="/home/somedir/user0/my_programs/TransDecoder-TransDecoder-v5.5.0";
    TDLONGORFS="${TDPATH}/TransDecoder.LongOrfs";
    TDPREDICT="${TDPATH}/TransDecoder.Predict";

    #Transdecoder Pfam, UniRef90, SwissProt DBs:
    TDPFAM="/disk0/user0/my_programs/transdecoder_dbs/Pfam-A.hmm";
    TDUR90="/disk0/user0/my_programs/transdecoder_dbs/transdecoder_dbs_diamond/uniref90.dmnd";
    #Using UniRef90 instead of SwissProt since it's more comprehensive.
    #
    TDSP="/disk0/user0/my_programs/transdecoder_dbs/transdecoder_dbs_diamond/uniprot_sprot.dmnd";


    #SAMPBNAME=$(basename ${SAMPFA} _trinity.Trinity.fasta);
    #SAMPBNAME=$(basename ${SAMPFA} _inf_mrna.fasta);
    SAMPBNAME=$(basename ${SAMPFA} _trinity.fasta);

    JOB_FILE="${TDJOBDIR}/${SAMPBNAME}_transdecoder_jobsub.bash";

    SAMPTDOUTDIR=${TDOUTDIR}/${SAMPBNAME};
    mkdir -p ${SAMPTDOUTDIR};


    #These are the outputs for the Blast and hmmer searches of the ORFs.
    BLASTOUT=${SAMPTDOUTDIR}/blastp.outfmt6;
    HMMOUT=${SAMPTDOUTDIR}/pfam.domtblout;

    #Longest ORFs file (need as input for homology and profile predictions).
    LORFS=${SAMPTDOUTDIR}/longest_orfs.pep;

    echo "#!/bin/bash -ex
#SBATCH -J minijob_td
#SBATCH -p hh
#SBATCH -o out.minijob.%J_${SAMPBNAME}
#SBATCH -e err.minijob.%J_${SAMPBNAME}
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=124G
#SBATCH --hint=nomultithread
#SBATCH -C "haswell"

mkdir -p ${SAMPTDOUTDIR};
cd ${SAMPTDOUTDIR};

#module load BIOPERL/1.6.923;

#Local copies of everything I need.
cp ${SAMPFA} ${SAMPTDOUTDIR}/$(basename ${SAMPFA});
cp ${TDSP} ${SAMPTDOUTDIR}/$(basename ${TDSP});
cp ${TDPFAM} ${SAMPTDOUTDIR}/$(basename ${TDPFAM});

#Longest ORFs prediction
#${TDLONGORFS} -t ${SAMPTDOUTDIR}/$(basename ${SAMPFA}) -O ${SAMPTDOUTDIR};
${TDLONGORFS} -t ${SAMPTDOUTDIR}/$(basename ${SAMPFA}) -O ${SAMPTDOUTDIR};


#Need this conda environment active since it contains
#the BLAST and hmmer installations that will be used by 
#TransDecoder.
#
source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
#conda activate transdecoder_conda;
#Homology predictions
#Blast predictions
#blastp -query ${LORFS} -db ${TDSP} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 -out ${BLASTOUT};
#conda deactivate;

#Using diamond instead of blast since blast is very slow.
#conda activate diamond_conda;
#diamond blastp --ultra-sensitive --query ${LORFS} --db ${TDUR90} --max-target-seqs 1 --outfmt 6 --evalue 1e-5 --threads 16 --out ${BLASTOUT};
#diamond blastp --ultra-sensitive --query ${LORFS} --db ${SAMPTDOUTDIR}/$(basename ${TDSP}) --max-target-seqs 1 --outfmt 6 --evalue 1e-5 --threads 16 --out ${BLASTOUT};
#conda deactivate;

#conda activate transdecoder_conda;
#HMM predictions
#hmmscan --cpu 16 --domtblout ${HMMOUT} ${TDPFAM} ${LORFS};
#hmmsearch --cpu 16 -o ${SAMPTDOUTDIR}/hmmsearch_outfile_irrelevant.tmp --domtblout ${SAMPTDOUTDIR}/hmmsearch.tmp ${SAMPTDOUTDIR}/$(basename ${TDPFAM}) ${LORFS};
#awk 'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=\$1; \$1=\$4; \$4=tmp; tmp=\$2; \$2=\$5; \$5=tmp; print}' ${SAMPTDOUTDIR}/hmmsearch.tmp > ${HMMOUT};
#conda deactivate;

#Transdecoder CDS predictions using these
#${TDPREDICT} --single_best_only -t ${SAMPTDOUTDIR}/$(basename ${SAMPFA}) -O ${SAMPTDOUTDIR} --retain_pfam_hits ${HMMOUT} --retain_blastp_hits ${BLASTOUT};
${TDPREDICT} --single_best_only -t ${SAMPTDOUTDIR}/$(basename ${SAMPFA}) -O ${SAMPTDOUTDIR};
#${TDPREDICT} -t ${SAMPFA} -O ${SAMPTDOUTDIR}


#Copying and renaming output files to make them more manageable.
cp ${SAMPTDOUTDIR}/$(basename ${SAMPFA}).transdecoder.cds ${MAINTDOUTDIR}/${SAMPBNAME}_transdecoder_cds.fasta;
cp ${SAMPTDOUTDIR}/$(basename ${SAMPFA}).transdecoder.pep ${MAINTDOUTDIR}/${SAMPBNAME}_transdecoder.fasta;


    " > ${JOB_FILE};
    sbatch ${JOB_FILE};

done

