#!/usr/bin/env bash
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=EDTA_
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/EDTA_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/EDTA_%J.err
#SBATCH --partition=pibu_el8


#Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/organization_annotation_course"
OUTDIR="$WORKDIR/results/EDTA_annotation"
LOGDIR="$WORKDIR/log"
APPTAINERPATH="/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif"
CDS="/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated"
#Assembly path
HIFIASM="/data/users/vmuller/assembly_annotation_course/results/Pacbio/05_assembly_Hifiasm/HiFiasm_Lu1_primary.fa"

#Create the directory for the error and output file if not present
mkdir -p "$LOGDIR"

#Create the directory output if not present
mkdir -p "$OUTDIR"

#move to the output dir so the result are writen in the correct directory
cd "$OUTDIR"

#Run EDTA 
apptainer exec --bind /data "$APPTAINERPATH" EDTA.pl  \
    --genome $HIFIASM \
    --species others \
    --step all \
    --sensitive 1 \
    --cds $CDS \
    --anno 1 \
    --threads $SLURM_CPUS_PER_TASK

