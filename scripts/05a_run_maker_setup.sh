#!/usr/bin/env bash
#SBATCH --job-name=maker_ctl_setup
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/maker_ctl_setup_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/maker_ctl_setup_%J.err


#Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/organization_annotation_course"
ANNODIR="$WORKDIR/results/annotation"
LOGDIR="$WORKDIR/log"

#input files 
ASSEMBLY_DIR="/data/users/${USER}/assembly_annotation_course"
GENOME="$ASSEMBLY_DIR/results/Pacbio/05_assembly_Hifiasm/HiFiasm_Lu1_primary.fa"
TRINITY="$ASSEMBLY_DIR/results/RNASeq/04_assembly_trinity/04_assembly_trinity.Trinity.fasta"

# EDTA outputs
EDTA_DIR="$WORKDIR/results/EDTA_annotation"
PREFIX="HiFiasm_Lu1_primary"
TELIB="$EDTA_DIR/${PREFIX}.fa.mod.EDTA.TElib.fa"

# Container
IMG="/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif"

# Evidence for the peptide
COURSE="/data/courses/assembly-annotation-course/CDS_annotation/data"
TAIR_PEP="$COURSE/TAIR10_pep_20110103_representative_gene_model"
UNIPROT_PEP="$COURSE/uniprot_viridiplantae_reviewed.fa"
PTREP="$COURSE/PTREP20"

#Creating the output and log directory if not present
mkdir -p "$ANNODIR" "$LOGDIR"
cd "$ANNODIR"

# Create control files
apptainer exec \
  --bind "$WORKDIR","/data/users/${USER}/assembly_annotation_course","/data/courses/assembly-annotation-course" \
  "$IMG" maker -CTL




# Configure maker_opts.ctl

# Setup temp directory
SCR="${SCRATCH:-$ANNODIR/tmp}"
mkdir -p "$SCR"

# Function to escape special characters for sed
esc() { printf '%s' "$1" | sed 's/[&/\]/\\&/g'; }

# Change the setting of the config files
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
