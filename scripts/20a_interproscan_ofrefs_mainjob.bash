#!/bin/bash -ex
#SBATCH -J Mainjob
#SBATCH -p hh
#SBATCH -o out.mainjob.%J
#SBATCH -e err.mainjob.%J
#SBATCH -t 12:00:00
#SBATCH -C "haswell"

#INTERPROSCAN="/disk0/user0/my_programs/interproscan-5.46-81.0/interproscan.sh";
#
INTERPROSCAN="/disk0/user0/my_programs/interproscan-5.51-85.0/interproscan.sh";

#Paths and stuff
MYPATH="/disk0/user0/newrun";

MAINJOBDIR="${MYPATH}/jobfiles";
IPJOBDIR="${MAINJOBDIR}/interproscan_ofrefs/jobsubs";
mkdir -p ${IPJOBDIR};

MAINOUTDIR="${MYPATH}/outputs";
IPOUTDIR="${MAINOUTDIR}/interproscan_ofrefs";
mkdir -p ${IPOUTDIR};


#OrthoFinder SOIs extracted using my custom R script.
OFDIR="${MAINOUTDIR}/orthofinder";
SOIDIR="${OFDIR}/rscript_of_getsois_ofrefs";



for SAMPFA in ${SOIDIR}/*.fasta;
do


    #Extracting sample basename.
    SAMPBNAME=$(basename ${SAMPFA} .fasta);

    #This variable sets the output path + output filename simultaneously
    #and will be passed as the input argument to the --output-file-base parameter
    IPSBN="${IPOUTDIR}/${SAMPBNAME}_interpro_res";
    IPSTMP="${IPOUTDIR}/${SAMPBNAME}_interpro_tmp";
    mkdir -p ${IPSTMP};

    #IPSOPTS="--verbose --cpu 16 --disable-precalc --goterms --pathways --iprlookup --formats tsv,svg,gff3,json,html --applications pfam,phobius,panther,mobidblite,cdd,gene3d --enable-tsv-residue-annot";
    IPSOPTS="--verbose --cpu 16 --disable-precalc --goterms --pathways --iprlookup --formats tsv,svg,gff3,html --applications pfam,panther,cdd,smart,gene3d,phobius,prositepatterns,prositeprofiles --enable-tsv-residue-annot";


    JOBFILE="${IPJOBDIR}/${SAMPBNAME}_job.bash";


    echo "#!/bin/bash -ex
#SBATCH -J minijob
#SBATCH -p hh
#SBATCH -o out.minijob_ips.%J
#SBATCH -e err.minijob_ips.%J
#SBATCH -t 12-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=124G
#SBATCH -C "haswell"

source /home/somedir/user0/miniconda3/etc/profile.d/conda.sh;

conda activate interproscan_conda;

cd ${IPOUTDIR};

${INTERPROSCAN} ${IPSOPTS} --tempdir ${IPSTMP} --output-file-base ${IPSBN} --input ${SAMPFA};

conda deactivate;

" > ${JOBFILE};

    sbatch ${JOBFILE};

done

