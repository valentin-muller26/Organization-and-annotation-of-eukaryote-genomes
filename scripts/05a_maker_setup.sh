#!/usr/bin/env bash
#SBATCH --job-name=maker_ctl_setup
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/maker_ctl_setup_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/maker_ctl_setup_%J.err

set -euo pipefail

# ========================================
# Configuration - ATTENTION AUX CHEMINS !
# ========================================
WORKDIR="/data/users/${USER}/organization_annotation_course"
ANNODIR="$WORKDIR/results/annotation"
LOGDIR="$WORKDIR/log"

# Input files - VÉRIFIE CES CHEMINS !
# Option 1: Si tes fichiers sont dans assembly_annotation_course
ASSEMBLY_DIR="/data/users/${USER}/assembly_annotation_course"
GENOME="$ASSEMBLY_DIR/results/Pacbio/05_assembly_Hifiasm/HiFiasm_Lu1_primary.fa"
TRINITY="$ASSEMBLY_DIR/results/RNASeq/04_assembly_trinity/04_assembly_trinity.Trinity.fasta"

# Option 2: Si tu as copié/lié les fichiers dans organization_annotation_course
# GENOME="$WORKDIR/results/Pacbio/05_assembly_Hifiasm/HiFiasm_Lu1_primary.fa"
# TRINITY="$WORKDIR/results/RNASeq/04_assembly_trinity/04_assembly_trinity.Trinity.fasta"

# EDTA outputs
EDTA_DIR="$WORKDIR/results/EDTA_annotation"
PREFIX="HiFiasm_Lu1_primary"
TELIB="$EDTA_DIR/${PREFIX}.fa.mod.EDTA.TElib.fa"

# Container
IMG="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"

# Evidence (course resources)
COURSE="/data/courses/assembly-annotation-course/CDS_annotation/data"
TAIR_PEP="$COURSE/TAIR10_pep_20110103_representative_gene_model"
UNIPROT_PEP="$COURSE/uniprot_viridiplantae_reviewed.fa"
PTREP="$COURSE/PTREP20"

# ========================================
# Setup
# ========================================
mkdir -p "$ANNODIR" "$LOGDIR"
cd "$ANNODIR"

echo "[INFO] Creating MAKER control files..."

# Create control files
apptainer exec \
  --bind "$WORKDIR","/data/users/${USER}/assembly_annotation_course","/data/courses/assembly-annotation-course" \
  "$IMG" maker -CTL

# ========================================
# Verify files exist - CRUCIAL !
# ========================================
echo "[INFO] Verifying input files..."

[ -s maker_opts.ctl ] || { echo "ERROR: maker_opts.ctl not created"; exit 2; }

echo "Checking genome: $GENOME"
[ -s "$GENOME" ] || { echo "ERROR: genome missing: $GENOME"; exit 3; }

echo "Checking EDTA TE library: $TELIB"
[ -s "$TELIB" ] || { echo "ERROR: EDTA TE library not found: $TELIB"; exit 4; }

echo "Checking TAIR proteins: $TAIR_PEP"
[ -s "$TAIR_PEP" ] || { echo "ERROR: TAIR protein set missing: $TAIR_PEP"; exit 5; }

echo "Checking UniProt proteins: $UNIPROT_PEP"
[ -s "$UNIPROT_PEP" ] || { echo "ERROR: UniProt plant set missing: $UNIPROT_PEP"; exit 6; }

echo "Checking PTREP20: $PTREP"
[ -s "$PTREP" ] || { echo "ERROR: PTREP20 missing: $PTREP"; exit 7; }

# Check for Trinity transcripts
if [ -s "$TRINITY" ]; then
  EST_VAL="$TRINITY"
  echo "[INFO] Using EST evidence: $TRINITY"
else
  EST_VAL=""
  echo "[WARN] No Trinity transcripts found — leaving est= blank"
  echo "[WARN] Searched at: $TRINITY"
fi

# ========================================
# Configure maker_opts.ctl
# ========================================
echo "[INFO] Configuring maker_opts.ctl..."

# Setup temp directory
SCR="${SCRATCH:-$ANNODIR/tmp}"
mkdir -p "$SCR"

# Function to escape special characters for sed
esc() { printf '%s' "$1" | sed 's/[&/\]/\\&/g'; }

# Patch maker_opts.ctl
sed -i \
  -e "s|^genome=.*|genome=$(esc "$GENOME")|" \
  -e "s|^est=.*|est=$(esc "$EST_VAL")|" \
  -e "s|^protein=.*|protein=$(esc "$TAIR_PEP"),$(esc "$UNIPROT_PEP")|" \
  -e "s|^model_org=.*|model_org=|" \
  -e "s|^rmlib=.*|rmlib=$(esc "$TELIB")|" \
  -e "s|^repeat_protein=.*|repeat_protein=$(esc "$PTREP")|" \
  -e "s|^augustus_species=.*|augustus_species=arabidopsis|" \
  -e "s|^est2genome=.*|est2genome=1|" \
  -e "s|^protein2genome=.*|protein2genome=1|" \
  -e "s|^cpus=.*|cpus=1|" \
  -e "s|^alt_splice=.*|alt_splice=1|" \
  -e "s|^TMP=.*|TMP=$(esc "$SCR")|" \
  maker_opts.ctl
