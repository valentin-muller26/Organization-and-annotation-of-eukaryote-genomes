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

# Path to assembly BUSCO results
GENOME_BUSCO="/data/users/vmuller/assembly_annotation_course/results/Pacbio/06_Busco_hifiasm"
TRANSCRIPTOME_BUSCO="/data/users/vmuller/assembly_annotation_course/results/RNASeq/05_Busco_trinity"

# Create the directory for the error and output file if not present
mkdir -p "$LOGDIR"

# Create the directory output if not present
mkdir -p "$OUTDIR"

# Change to output directory
cd "$OUTDIR"

# Load BUSCO module
module load BUSCO/5.4.2-foss-2021a

# Generate individual plots
generate_plot.py -wd "$OUTDIR/busco_protein"
generate_plot.py -wd "$OUTDIR/busco_transcript"

# Combined plot with all results
mkdir -p "$OUTDIR/combined_summaries"

# Copy annotation files
cp "$OUTDIR/busco_protein"/short_summary*.txt "$OUTDIR/combined_summaries/" 
cp "$OUTDIR/busco_transcript"/short_summary*.txt "$OUTDIR/combined_summaries/" 

# Copy assembly files (searching in subdirectories)
find "$GENOME_BUSCO" -name "short_summary*.txt" -type f -exec cp {} "$OUTDIR/combined_summaries/" 
find "$TRANSCRIPTOME_BUSCO" -name "short_summary*.txt" -type f -exec cp {} "$OUTDIR/combined_summaries/" 

# Rename files in combined_summaries
cd "$OUTDIR/combined_summaries"

#Rename the file to have a good name in the horizontal track of the plot
# Rename busco_protein -> annotation_protein
for file in *busco_protein*; do
    [ -f "$file" ] && mv "$file" "${file//busco_protein/annotation_protein}"
done

# Rename busco_transcript -> annotation_transcript
for file in *busco_transcript*; do
    [ -f "$file" ] && mv "$file" "${file//busco_transcript/annotation_transcript}"
done

# Rename 06_Busco_hifiasm -> assembly_genome
for file in *06_Busco_hifiasm*; do
    [ -f "$file" ] && mv "$file" "${file//06_Busco_hifiasm/assembly_genome}"
done

# Rename 05_Busco_trinity -> assembly_transcriptome
for file in *05_Busco_trinity*; do
    [ -f "$file" ] && mv "$file" "${file//05_Busco_trinity/assembly_transcriptome}"
done

echo "=== Files after renaming ==="
ls -la

# Generate combined plot
cd "$OUTDIR"
generate_plot.py -wd "$OUTDIR/combined_summaries"

