#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=TEsorter
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/TEsorter_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/TEsorter_%J.err
#SBATCH --partition=pibu_el8

# ============================================================================
# STEP 2: CLASSIFY LTR-RTs WITH TESORTER
# ============================================================================

# Directories
WORKDIR="/data/users/${USER}/organization_annotation_course"
EDTA_DIR="$WORKDIR/results/EDTA_annotation"
OUTDIR="$WORKDIR/results/LTR_analysis"

# Container
TESORTER_CONTAINER="/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif"

# Input file
LTR_FASTA="$EDTA_DIR/HiFiasm_Lu1_primary.fa.mod.EDTA.raw/LTR/HiFiasm_Lu1_primary.fa.mod.LTR.raw.fa"

# Create output directory
mkdir -p "$OUTDIR"

# Change to output directory
cd "$OUTDIR"

# Run TEsorter
apptainer exec \
    --bind "$WORKDIR" \
    "$TESORTER_CONTAINER" \
    TEsorter \
        "$LTR_FASTA" \
        -db rexdb-plant \
        -p "$SLURM_CPUS_PER_TASK"

