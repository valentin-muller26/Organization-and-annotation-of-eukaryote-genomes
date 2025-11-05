#!/usr/bin/env bash
#SBATCH --job-name=genespace
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/genespace_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/genespace_%J.err

# Directories
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/vmuller/organization_annotation_course"
GENESPACE_WD="${WORKDIR}/results/genespace"

# Container and script
CONTAINER="${COURSEDIR}/containers/genespace_latest.sif"
RSCRIPT="${WORKDIR}/scripts/08b_genespace.R"

# Run GENESPACE in container
echo "Starting GENESPACE analysis..."
echo ""

apptainer exec \
    --bind /data \
    --bind ${SCRATCH}:/temp \
    ${CONTAINER} Rscript ${RSCRIPT} ${GENESPACE_WD}

echo ""
echo "=========================================="
echo "Job completed successfully!"
echo "=========================================="