#!/usr/bin/env bash
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_index
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/fai_index_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/fai_index_%J.err
#SBATCH --partition=pibu_el8


#Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/organization_annotation_course"
OUTDIR="$WORKDIR/results/fai_files"
LOGDIR="$WORKDIR/log"

#Assembly path
HIFIASM="/data/users/vmuller/assembly_annotation_course/results/Pacbio/05_assembly_Hifiasm/HiFiasm_Lu1_primary.fa"

#Create the directory for the error and output file if not present
mkdir -p "$LOGDIR"

#Create the directory output if not present
mkdir -p "$OUTDIR"

#Move to the output dir so the result are written in the correct directory
cd "$OUTDIR"

#Load samtools module
module load SAMtools/1.13-GCC-10.3.0

#Generate FAI index file
samtools faidx "$HIFIASM"

#Copy the FAI file to the output directory
cp "${HIFIASM}.fai" "$OUTDIR/"
