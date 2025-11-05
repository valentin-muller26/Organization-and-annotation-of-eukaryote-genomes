#!/usr/bin/env bash
#SBATCH --job-name=Busco
#SBATCH --partition=pibu_el8
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/Busco_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/Busco_%J.err

WORKDIR="/data/users/${USER}/organization_annotation_course"
LOGDIR="$WORKDIR/log"

# Path to the directory of the input files
ANNODIR="$WORKDIR/results/annotation"
INPUTDIR="$ANNODIR/final"
OUTDIR="$ANNODIR/busco"

#input fasta files
PROTEINFILE="HiFiasm_Lu1_primary.all.maker.proteins.renamed.filtered.longest2.fasta"
TRANSCRIPTFILE="HiFiasm_Lu1_primary.all.maker.transcripts.renamed.filtered.longest2.fasta"

# Create the directory for the error and output file if not present
mkdir -p "$LOGDIR"

# Create the directory output if not present
mkdir -p "$OUTDIR"

# Change to output directory
cd "$OUTDIR"

# Load BUSCO module
module load BUSCO/5.4.2-foss-2021a

# Run BUSCO on proteins
busco -i "$INPUTDIR/$PROTEINFILE" \
      -l brassicales_odb10 \
      -o busco_protein \
      -m proteins \
      -c 4


# Run BUSCO on transcripts
busco -i "$INPUTDIR/$TRANSCRIPTFILE" \
      -l brassicales_odb10 \
      -o busco_transcript \
      -m transcriptome \
      -c 4

