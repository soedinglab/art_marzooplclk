#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00
#SBATCH -C "haswell"


#Paths and stuff
MYPATH="/disk0/user0/newrun";
#EMAPPER="/disk0/user0/my_programs/eggnog-mapper/emapper.py";
EMAPPER="/disk0/user0/my_programs/eggnog-mapper-2.1.4-2/emapper.py";

MAINJOBDIR="${MYPATH}/jobfiles";
EMAPJOBDIR="${MAINJOBDIR}/eggnogmapper/jobsubs";
mkdir -p ${EMAPJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
EMAPOUTDIR="${MAINOUTDIR}/eggnogmapper";
mkdir -p ${EMAPOUTDIR};



#Path to eggnogmapper databases directory.
#Has Diamond and Hmmer + Pfam databases.
EMAPDBPATH="/disk0/user0/my_programs/eggnog-mapper/dbs";

#Assemblies path
#
ASSEMPATH="${MAINOUTDIR}/rscript_lenfilt";
#ASSEMPATH="${MAINOUTDIR}/mmseqs2_protclu_trinity";

#
#for SAMPFA in ${ASSEMPATH}/*_pclu_rep_seq.fasta;
for SAMPFA in ${ASSEMPATH}/*_td_rfilt.fasta; 
do

    #
    SAMPBNAME=$(basename ${SAMPFA} _td_rfilt.fasta);
    #SAMPBNAME=$(basename ${SAMPFA} _pclu_rep_seq.fasta);
    JOB_FILE="${EMAPJOBDIR}/${SAMPBNAME}_eggnogmapper_job.bash";

    SAMPOUTDIR="${EMAPOUTDIR}/${SAMPBNAME}_eggnogmapper";


    echo "#!/bin/bash -ex
#SBATCH -J minijob_iso
#SBATCH -p hh
#SBATCH -o out.minijob_iso.%J
#SBATCH -e err.minijob_iso.%J
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=124G
#SBATCH -C "haswell"

#
source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;
conda activate eggnogmapperinstall_conda;

#export PATH=/disk0/user0/my_programs/eggnog-mapper-2.1.2/:/disk0/user0/my_programs/eggnog-mapper-2.1.2/eggnogmapper/bin:"\$PATH";
export PATH=/disk0/user0/my_programs/eggnog-mapper-2.1.4-2/:/disk0/user0/my_programs/eggnog-mapper-2.1.4-2/eggnogmapper/bin:"\$PATH";

mkdir -p ${SAMPOUTDIR};
cd  ${SAMPOUTDIR};

#emapper.py
${EMAPPER} --output_dir ${SAMPOUTDIR} --output ${SAMPBNAME} --pfam_realign realign --go_evidence all --report_orthologs --tax_scope Metazoa --sensmode ultra-sensitive -m diamond --itype proteins --cpu 16 --data_dir ${EMAPDBPATH} -i ${SAMPFA};

conda deactivate;

mv ${SAMPOUTDIR}/* ${EMAPOUTDIR};

    " > ${JOB_FILE};
   sbatch ${JOB_FILE};

done

