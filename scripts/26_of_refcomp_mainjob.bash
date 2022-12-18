#!/bin/bash -ex
#SBATCH -J mainjob_special
#SBATCH -p hh
#SBATCH -o out.mainjob_special.%J
#SBATCH -e err.mainjob_special.%J
#SBATCH -t 14-00:00:00
#SBATCH --cpus-per-task=256
#SBATCH --mem=0

#This is a bit of a manual script since it's executing just one single giant assembly that's executed off of 
#the compute node's local NVMe drive.

#Paths and stuff
MYPATH="/disk0/user0/newrun";

#Jobs
MAINJOBDIR="${MYPATH}/jobfiles";
SPJOBDIR="${MAINJOBDIR}/of_refcomp/jobsubs";
mkdir -p ${SPJOBDIR};

#Outputs
MAINOUTDIR="${MYPATH}/outputs";
MAINSPOUTDIR="${MAINOUTDIR}/of_refcomp";
mkdir -p ${MAINSPOUTDIR};

NVMEPATH="/nvme/n00/user0/newrun";
SPOUTDIR="${NVMEPATH}/outputs/of_refcomp";
mkdir -p ${SPOUTDIR};

#OrthoFinder needs the translated assemblies and the reference proteomes
#carrying the sequences of interest to which I wish to find orthologs.
ASSEMPATH="${MAINOUTDIR}/rscript_lenfilt";
REFASSEM="${MYPATH}/cc_orthofinder_dbs/final_ref_orthofinder";
REFASSEM_EXT="${MYPATH}/cc_orthofinder_dbs/extra_metazoan_refs";


#I will copy all of these to ${SPOUTDIR}/inputs.
cd ${SPOUTDIR};
mkdir -p ${SPOUTDIR}/inputs;
#
echo "\n\nCopying over TARGETS to inputs directory and renaming!!\n\n";
cp ${ASSEMPATH}/*_td_rfilt.fasta ${SPOUTDIR}/inputs;

echo "\n\nCopying over REFERENCES to inputs directory and renaming!!\n\n";
#
cp ${REFASSEM}/*_allrefs.fasta ${SPOUTDIR}/inputs;

echo "\n\nCopying over EXTRA REFERENCES to inputs directory and renaming!!\n\n";
cp ${REFASSEM_EXT}/*.fasta ${SPOUTDIR}/inputs;

echo "\n\nEXECUTING ORTHOFINDER!!\n\n";

#OrthoFinder works better when the user supplies their own species tree.
#I have one made from NCBI taxid data.
OFRSPTREE="${MYPATH}/cc_orthofinder_dbs/of_sp_tree_extrefs.nwk";
LOCSPTREE="${SPOUTDIR}/of_sp_tree_extrefs.nwk";
cp ${OFRSPTREE} ${LOCSPTREE};


#OrthoFinder run.
source /home/somedir/user0/my_programs/miniconda3/etc/profile.d/conda.sh;
conda activate orthofinder_conda;

orthofinder -t 512 -a 512 -M msa -S diamond_ultra_sens -s ${LOCSPTREE} -y -o ${SPOUTDIR}/outputs -f ${SPOUTDIR}/inputs;

conda deactivate;

#Copying over the entire OrthoFinder directory.
echo "\n\nCOPYING OVER OUTPUT TO DESTINATION DIRECTORY!!\n\n";
cp -r ${SPOUTDIR}/* ${MAINSPOUTDIR};

echo "\n\nALL DONE!!\n\n";
