#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=120G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=Maker
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/Maker_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/Maker_%J.err

set -euo pipefail

# ========================================
# Configuration
# ========================================
WORKDIR="/data/users/${USER}/organization_annotation_course"
ASSEMBLY_DIR="/data/users/${USER}/assembly_annotation_course"
ANNODIR="$WORKDIR/results/annotation"
LOGDIR="$WORKDIR/log"

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
REPEATMASKER_DIR="$COURSEDIR/softwares/RepeatMasker"
IMG="$COURSEDIR/containers/MAKER_3.01.03.sif"

export PATH="$PATH:$REPEATMASKER_DIR"

# ========================================
# Load modules
# ========================================
module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

# ========================================
# Verify setup
# ========================================
mkdir -p "$LOGDIR"
cd "$ANNODIR" || { echo "ERROR: cannot cd to $ANNODIR"; exit 1; }

echo "[INFO] Working directory: $ANNODIR"
echo "[INFO] Using $SLURM_NTASKS_PER_NODE MPI tasks"

# Check control files exist
for f in maker_opts.ctl maker_bopts.ctl maker_exe.ctl; do
  [ -s "$f" ] || { echo "ERROR: $f not found or empty"; exit 2; }
done

# ========================================
# Run MAKER with MPI
# ========================================
echo "[INFO] Starting MAKER annotation..."
echo "[INFO] Start time: $(date)"

mpiexec --oversubscribe -n $SLURM_NTASKS_PER_NODE apptainer exec \
  --bind /data \
  --bind "$SCRATCH":/TMP \
  --bind "$COURSEDIR" \
  --bind "$AUGUSTUS_CONFIG_PATH" \
  --bind "$REPEATMASKER_DIR" \
  "$IMG" \
  maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_exe.ctl

# ========================================
# Post-run summary
# ========================================
echo
echo "✅ MAKER run finished"
echo "End time: $(date)"
echo "Output directory: $ANNODIR"
echo "Logs: $LOGDIR"