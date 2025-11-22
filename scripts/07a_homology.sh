#!/usr/bin/env bash
#SBATCH --job-name=homology
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --output=/data/users/vmuller/organization_annotation_course/log/homology_%J.out
#SBATCH --error=/data/users/vmuller/organization_annotation_course/log/homology_%J.err

# General path
WORKDIR="/data/users/${USER}/organization_annotation_course"
LOGDIR="$WORKDIR/log"
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

# Path to the directory of the input files
ANNODIR="$WORKDIR/results/annotation/final"
FASTAPROTEINFILE="HiFiasm_Lu1_primary.all.maker.proteins.fasta.renamed.filtered.fasta"
GFFFILE="filtered.genes.renamed.gff3"
OUTDIR="$WORKDIR/results/homology"
mkdir -p $OUTDIR
cd $OUTDIR

# Path to the program used to merge
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"

# Define databases
UNIPROTDB="$COURSEDIR/data/uniprot/uniprot_viridiplantae_reviewed.fa"
TAIR10DB="$COURSEDIR/data/TAIR10_pep_20110103_representative_gene_model"

# Load module
module load BLAST+/2.15.0-gompi-2021a

################################################################################
# BLAST against UniProt Viridiplantae
################################################################################
echo "Running BLAST against UniProt..."
BLASTP_UNIPROT="blastp_uniprot.out"

blastp -query $ANNODIR/$FASTAPROTEINFILE \
    -db $UNIPROTDB \
    -num_threads 12 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out $BLASTP_UNIPROT

# Sort to keep only the best hit per query sequence
echo "Sorting UniProt BLAST results..."
sort -k1,1 -k12,12gr $BLASTP_UNIPROT | sort -u -k1,1 --merge > ${BLASTP_UNIPROT}.besthits

# Copy original files before adding functional annotations
echo "Creating backup copies..."
cp $ANNODIR/$FASTAPROTEINFILE ${FASTAPROTEINFILE}.Uniprot
cp $ANNODIR/$GFFFILE ${GFFFILE}.Uniprot

# Map the protein putative functions to the MAKER produced GFF3 and FASTA files
echo "Mapping UniProt functions to MAKER files..."
$MAKERBIN/maker_functional_fasta $UNIPROTDB ${BLASTP_UNIPROT}.besthits \
    ${FASTAPROTEINFILE}.Uniprot > ${FASTAPROTEINFILE}.Uniprot.annotated

$MAKERBIN/maker_functional_gff $UNIPROTDB ${BLASTP_UNIPROT}.besthits \
    ${GFFFILE}.Uniprot > ${GFFFILE}.Uniprot.annotated

################################################################################
# BLAST against Arabidopsis thaliana TAIR10
################################################################################
echo "Running BLAST against TAIR10..."
BLASTP_TAIR10="blastp_tair10.out"

blastp -query $ANNODIR/$FASTAPROTEINFILE \
    -db $TAIR10DB \
    -num_threads 12 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out $BLASTP_TAIR10

# Sort to keep only the best hit per query sequence
echo "Sorting TAIR10 BLAST results..."
sort -k1,1 -k12,12gr $BLASTP_TAIR10 | sort -u -k1,1 --merge > ${BLASTP_TAIR10}.besthits
