#!/bin/bash
#SBATCH --job-name=extract_TE_families
#SBATCH --partition=pibu_el8
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/extract_TE_families_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/extract_TE_families_%J.err

# Defining the constant for the path and files
WORKDIR="/data/users/${USER}/organization_annotation_course/results/EDTA_annotation"
GENOME="HiFiasm_Lu1_primary.fa"
TELIB="${WORKDIR}/${GENOME}.mod.EDTA.TElib.fa"

#Load the module for SeqKit 
module load SeqKit/2.6.1

cd "$WORKDIR"

# Extract Copia sequences
seqkit grep -r -p "Copia" "$TELIB" > Copia_sequences.fa

# Extract Gypsy sequences
seqkit grep -r -p "Gypsy" "$TELIB" > Gypsy_sequences.fa
