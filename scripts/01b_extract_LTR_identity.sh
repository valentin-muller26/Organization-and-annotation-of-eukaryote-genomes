#!/usr/bin/env bash
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Extract_LTR_identity
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/extract_identity_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/extract_identity_%J.err
#SBATCH --partition=pibu_el8

# ============================================================================
# STEP 1: EXTRACT LTR PERCENT IDENTITY FROM GFF3
# ============================================================================

# Directories
WORKDIR="/data/users/${USER}/organization_annotation_course"
EDTA_DIR="$WORKDIR/results/EDTA_annotation"
OUTDIR="$WORKDIR/results/LTR_analysis"

# Input file
GFF3_FILE="$EDTA_DIR/HiFiasm_Lu1_primary.fa.mod.EDTA.raw/LTR/HiFiasm_Lu1_primary.fa.mod.LTR.intact.raw.gff3"

# Create output directory
mkdir -p "$OUTDIR"

# Extract LTR identity from GFF3
# The GFF3 contains attributes like "ltr_identity=98.5"
grep -v "^#" "$GFF3_FILE" | \
    awk -F'\t' '
    /LTR_retrotransposon/ {
        # Extract the ID and ltr_identity from attributes
        match($9, /ID=([^;]+)/, id_arr)
        match($9, /ltr_identity=([0-9.]+)/, identity_arr)
        
        if (length(identity_arr) > 0) {
            print id_arr[1] "\t" identity_arr[1]
        }
    }' > "$OUTDIR/LTR_identity_raw.tsv"

# Count extracted entries
n_entries=$(wc -l < "$OUTDIR/LTR_identity_raw.tsv")
echo "Extracted $n_entries LTR-RTs with identity information"

if [ "$n_entries" -eq 0 ]; then
    echo "WARNING: No LTR identity data found in GFF3"
    echo "Trying alternative extraction method..."
    
    # Alternative method: look for any numeric value after ltr_identity
    grep "ltr_identity" "$GFF3_FILE" | \
        grep -oP 'ltr_identity=[0-9.]+' | \
        sed 's/ltr_identity=//' | \
        awk '{print "LTR_" NR "\t" $1}' > "$OUTDIR/LTR_identity_raw.tsv"
    
    n_entries=$(wc -l < "$OUTDIR/LTR_identity_raw.tsv")
    echo "Extracted $n_entries entries using alternative method"
fi

# Calculate basic statistics
echo ""
echo "LTR Identity Statistics:"
awk '{
    sum+=$2; 
    count++; 
    if(min==""){min=max=$2}; 
    if($2>max){max=$2}; 
    if($2<min){min=$2}
} 
END {
    print "  Count: " count
    print "  Mean: " sum/count "%"
    print "  Min: " min "%"
    print "  Max: " max "%"
}' "$OUTDIR/LTR_identity_raw.tsv"

