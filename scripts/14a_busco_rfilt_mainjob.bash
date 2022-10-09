#!/bin/bash -ex
#SBATCH -J Rerun
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 96:00:00
#SBATCH --hint=nomultithread
#SBATCH -C "haswell"

#module load JAVA/OPENJDK/13.0.1


#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
BUSCOJOBDIR="${MAINJOBDIR}/busco_rfilt/jobsubs";
mkdir -p ${BUSCOJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
BUSCOOUTDIR="${MAINOUTDIR}/busco_rfilt";
mkdir -p ${BUSCOOUTDIR};


#Assemblies path -- these are the transcriptomes filtered
#on the basis of rfilt proteomes.
ASSEMPATH=${MAINOUTDIR}//rscript_annotscmp/transcriptomes;


for SAMPFA in ${ASSEMPATH}/*_transcriptome.fasta;
do
    
    SAMPBNAME=$(basename ${SAMPFA} _transcriptome.fasta);

    #File name shenanigans since the rfilt transcriptomes do not have the grp suffix in them.
    if [ "${SAMPBNAME}" == "Acartia_tonsa" ] || [ "${SAMPBNAME}" == "Acartia_clausii" ] || [ "${SAMPBNAME}" == "Evadne_nordmanni" ] || [ "${SAMPBNAME}" == "Oikopleura_dioica" ] || [ "${SAMPBNAME}" == "Temora_longicornis" ];
    then
        SAMPBNAME="${SAMPBNAME}grp";
    fi;

    JOB_FILE="${BUSCOJOBDIR}/${SAMPBNAME}_busco_trinity_jobsub.bash";

    #Busco databases
    BUSCODB="/disk0/user0/my_programs/busco_downloads/lineages";

    

echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob.%J
#SBATCH -e err.minijob.%J
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=7G
#SBATCH --hint=nomultithread
#SBATCH --exclude=hh[002-003]
#SBATCH -C "haswell"

source /usr/users/user0/miniconda3/etc/profile.d/conda.sh
conda activate busco_conda

cd ${BUSCOOUTDIR};

#Eukaryota
busco --in ${SAMPFA} --out ${SAMPBNAME}_busco_euk --out_path ${BUSCOOUTDIR} --mode transcriptome --cpu 16 --offline --lineage ${BUSCODB}/eukaryota_odb10;

#Arthropoda since most samples are arthropods
busco --in ${SAMPFA} --out ${SAMPBNAME}_busco_arth --out_path ${BUSCOOUTDIR} --mode transcriptome --cpu 16 --offline --lineage ${BUSCODB}/arthropoda_odb10;


conda deactivate

    " > ${JOB_FILE};
sbatch ${JOB_FILE};
done



