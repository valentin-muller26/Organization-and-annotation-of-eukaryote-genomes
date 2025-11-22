#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=120G
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50
#SBATCH --job-name=Maker
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/Maker_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/Maker_%J.err


#Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/organization_annotation_course"
ASSEMBLY_DIR="/data/users/${USER}/assembly_annotation_course"
ANNODIR="$WORKDIR/results/annotation"
LOGDIR="$WORKDIR/log"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
REPEATMASKER_DIR="$COURSEDIR/softwares/RepeatMasker"
#Container 
IMG="$COURSEDIR/containers/MAKER_3.01.03.sif"

export PATH="$PATH:$REPEATMASKER_DIR"


# Load modules
module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

#Creating log dir if not exist and moving to working directory
mkdir -p "$LOGDIR"
cd "$ANNODIR" 

# Run MAKER with MPI
mpiexec --oversubscribe -n $SLURM_NTASKS_PER_NODE apptainer exec \
  --bind /data \
  --bind "$SCRATCH":/TMP \
  --bind "$COURSEDIR" \
  --bind "$AUGUSTUS_CONFIG_PATH" \
  --bind "$REPEATMASKER_DIR" \
  "$IMG" \
  maker -mpi --ignore_nfs_tmp -TMP /TMP maker_opts.ctl maker_bopts.ctl maker_exe.ctl
