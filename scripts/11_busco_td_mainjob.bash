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
BUSCOJOBDIR="${MAINJOBDIR}/busco_td_trinity/jobsubs";
mkdir -p ${BUSCOJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
BUSCOOUTDIR="${MAINOUTDIR}/busco_td_trinity";
mkdir -p ${BUSCOOUTDIR};


#Assemblies path
ASSEMPATH=${MAINOUTDIR}/transdecoder_trinity;


for SAMPFA in ${ASSEMPATH}/*_transdecoder_fixed_noast.fasta;
do
    
    SAMPBNAME=$(basename ${SAMPFA} _transdecoder_fixed_noast.fasta);

    JOB_FILE="${BUSCOJOBDIR}/${SAMPBNAME}_busco_td_trinity_jobsub.bash";

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

source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh
conda activate busco_conda

cd ${BUSCOOUTDIR};

#Eukaryota
busco --in ${SAMPFA} --out ${SAMPBNAME}_busco_euk_tdtrinity --out_path ${BUSCOOUTDIR} --mode proteins --cpu 16 --offline --lineage ${BUSCODB}/eukaryota_odb10;

#Arthropoda since most samples are arthropods
busco --in ${SAMPFA} --out ${SAMPBNAME}_busco_arth_tdtrinity --out_path ${BUSCOOUTDIR} --mode proteins --cpu 16 --offline --lineage ${BUSCODB}/arthropoda_odb10;


conda deactivate

    " > ${JOB_FILE};
sbatch ${JOB_FILE};
done



