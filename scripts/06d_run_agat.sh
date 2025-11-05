#!/usr/bin/env bash
#SBATCH --job-name=AGAT_stats
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/AGAT_stats_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/AGAT_stats_%J.err

WORKDIR="/data/users/${USER}/organization_annotation_course"
LOGDIR="$WORKDIR/log"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

# Path to the directory of the input files
ANNODIR="$WORKDIR/results/annotation/"
INPUTDIR="$ANNODIR/final"
INPUTFILE="filtered.genes.renamed.gff3"
OUTDIR="$ANNODIR/agat_stats"

APPTAINERPATH="/data/courses/assembly-annotation-course/CDS_annotation/containers/agat_1.5.1--pl5321hdfd78af_0.sif"

# Create directories if not present
mkdir -p "$LOGDIR"
mkdir -p "$OUTDIR"

# Run AGAT statistics using apptainer
apptainer exec --bind /data "$APPTAINERPATH" agat_sp_statistics.pl \
  -i "$INPUTDIR/$INPUTFILE" \
  -o "$OUTDIR/annotation_statistics.txt"