#!/bin/bash
#SBATCH --job-name=TEsorter_Copia_Gypsy
#SBATCH --partition=pibu_el8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=10
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/TEsorter_Copia_Gypsy_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/TEsorter_Copia_Gypsy_%J.err

# Configuration
WORKDIR="/data/users/${USER}/organization_annotation_course"
INPUTDIR="$WORKDIR/results/EDTA_annotation"
CONTAINER="/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif"
OUTDIR="$WORKDIR/results/TEsorter_Copia_Gypsy"
LOGDIR="$WORKDIR/log"

# Create directories
mkdir -p "$LOGDIR" "$OUTDIR"

# Copia classification
apptainer exec --bind /data "$CONTAINER" \
  TEsorter $INPUTDIR/Copia_sequences.fa -db rexdb-plant -p "$SLURM_CPUS_PER_TASK" -pre "$OUTDIR/Copia"

# Gypsy classification
apptainer exec --bind /data "$CONTAINER" \
  TEsorter $INPUTDIR/Gypsy_sequences.fa -db rexdb-plant -p "$SLURM_CPUS_PER_TASK" -pre "$OUTDIR/Gypsy"