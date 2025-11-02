#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8
#SBATCH --job-name=extract_longest2
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/extract_longest2_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/extract_longest2_%J.err


# General path
WORKDIR="/data/users/${USER}/organization_annotation_course"
ANNODIR="$WORKDIR/results/annotation/final"

# Change to final directory
cd "$ANNODIR" || exit 1

# Load samtools
module load SAMtools/1.13-GCC-10.3.0

# Define file names
protein="HiFiasm_Lu1_primary.all.maker.proteins.fasta.renamed.filtered.fasta"
transcript="HiFiasm_Lu1_primary.all.maker.transcripts.fasta.renamed.filtered.fasta"

echo "========================================"
echo "Input files:"
echo "  Protein: $protein"
echo "  Transcript: $transcript"
echo "========================================"

# Check if input files exist
if [ ! -f "$protein" ]; then
    echo "ERROR: Protein file not found: $protein"
    exit 1
fi

if [ ! -f "$transcript" ]; then
    echo "ERROR: Transcript file not found: $transcript"
    exit 1
fi

# Extract Longest Protein Isoforms using samtools
echo ""
echo "Extracting longest protein isoforms with samtools..."

# Step 1: Index the fasta file
samtools faidx "$protein"

# Step 2: Get sequence lengths and extract gene names
cut -f1,2 ${protein}.fai | \
awk '{
    # Extract gene name (everything before -R)
    gene = $1
    sub(/-R.*/, "", gene)
    # Keep only longest isoform per gene
    if ($2 > maxlen[gene]) {
        maxlen[gene] = $2
        longest[gene] = $1
    }
}
END {
    for (gene in longest) {
        print longest[gene]
    }
}' > protein_longest_ids.txt

# Step 3: Extract sequences using samtools
samtools faidx "$protein" $(cat protein_longest_ids.txt | tr '\n' ' ') > HiFiasm_Lu1_primary.all.maker.proteins.renamed.filtered.longest2.fasta

# Extract Longest Transcript Isoforms using samtools
echo "Extracting longest transcript isoforms with samtools..."

# Step 1: Index the fasta file
samtools faidx "$transcript"

# Step 2: Get sequence lengths and extract gene names
cut -f1,2 ${transcript}.fai | \
awk '{
    # Extract gene name (everything before -R)
    gene = $1
    sub(/-R.*/, "", gene)
    # Keep only longest isoform per gene
    if ($2 > maxlen[gene]) {
        maxlen[gene] = $2
        longest[gene] = $1
    }
}
END {
    for (gene in longest) {
        print longest[gene]
    }
}' > transcript_longest_ids.txt

# Step 3: Extract sequences using samtools
samtools faidx "$transcript" $(cat transcript_longest_ids.txt | tr '\n' ' ') > HiFiasm_Lu1_primary.all.maker.transcripts.renamed.filtered.longest2.fasta

# Verify counts and file sizes
echo ""
echo "========================================"
echo "Results:"
echo "========================================"
echo "Original protein sequences:    $(grep -c '^>' $protein)"
echo "Longest protein isoforms:      $(grep -c '^>' HiFiasm_Lu1_primary.all.maker.proteins.renamed.filtered.longest2.fasta)"
echo "Protein output size:           $(du -h HiFiasm_Lu1_primary.all.maker.proteins.renamed.filtered.longest2.fasta | cut -f1)"
echo ""
echo "Original transcript sequences: $(grep -c '^>' $transcript)"
echo "Longest transcript isoforms:   $(grep -c '^>' HiFiasm_Lu1_primary.all.maker.transcripts.renamed.filtered.longest2.fasta)"
echo "Transcript output size:        $(du -h HiFiasm_Lu1_primary.all.maker.transcripts.renamed.filtered.longest2.fasta | cut -f1)"
echo ""
echo "Files created:"
echo "  - HiFiasm_Lu1_primary.all.maker.proteins.renamed.filtered.longest2.fasta"
echo "  - HiFiasm_Lu1_primary.all.maker.transcripts.renamed.filtered.longest2.fasta"
echo ""
echo "Extraction complete!"

# Clean up intermediate files
rm -f protein_longest_ids.txt transcript_longest_ids.txt
rm -f ${protein}.fai ${transcript}.fai