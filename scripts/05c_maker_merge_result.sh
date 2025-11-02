#!/usr/bin/env bash
#SBATCH --job-name=maker_mergeresult
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/maker_mergeresult_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/maker_mergeresult_%J.err

#General  path
WORKDIR="/data/users/${USER}/organization_annotation_course"
LOGDIR="$WORKDIR/log"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

#path to the directory of the input files
ANNODIR="$WORKDIR/results/annotation/"
DATASTORE="$ANNODIR/HiFiasm_Lu1_primary.maker.output"
MASTERINDEXFILE="HiFiasm_Lu1_primary_master_datastore_index.log"

#Path to the program used to merge
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"
# Create output directory if it doesn't exist
mkdir -p "$ANNODIR/Merge_result"

# Merge GFF with sequences
$MAKERBIN/gff3_merge -s -d $DATASTORE/$MASTERINDEXFILE > "$ANNODIR/Merge_result/HiFiasm_Lu1_primary.all.maker.gff"

# Merge GFF without sequences
$MAKERBIN/gff3_merge -n -s -d $DATASTORE/$MASTERINDEXFILE > "$ANNODIR/Merge_result/HiFiasm_Lu1_primary.all.maker.noseq.gff"

# Merge FASTA files
$MAKERBIN/fasta_merge -d $DATASTORE/$MASTERINDEXFILE -o "$ANNODIR/Merge_result/HiFiasm_Lu1_primary"
