# Script: MGFseq_mapping-counting.sh
# Description:
# This script performs read mapping, BAM file processing, and feature counting for MGFseq data.
# It uses minimap2 for mapping, samtools for BAM processing, and featureCounts for counting mapped reads.
# Author: Lissa Cruz-Saavedra
# Date: 02-04-2025
# Version: 1.0

# Exit immediately if a command exits with a non-zero status
set -e

# Variables (update these as needed)
DATA="data"                          # Directory containing input FASTQ files
GENOME="genome.fasta"                # Reference genome file
GTF="annotations.gtf"                # Annotation file for featureCounts
THREADS=4                            # Number of threads to use
SAMPLEDIR="sample_output"            # Output directory for results

# Create necessary directories
mkdir -p "${SAMPLEDIR}/result_mapcount/sam"
mkdir -p "${SAMPLEDIR}/result_mapcount/bam"
mkdir -p "${SAMPLEDIR}/result_mapcount/sorted_bam"
mkdir -p "${SAMPLEDIR}/result_mapcount/counts"

# Step 1: Mapping reads using minimap2
echo "Step 1: Mapping reads using minimap2..."
for sample in "$DATA"/*.fastq; do
    sample_name=$(basename "$sample" .fastq)
    echo "Mapping $sample_name..."
    minimap2 -ax map-ont -t "$THREADS" "$GENOME" "$sample" > "${SAMPLEDIR}/result_mapcount/sam/${sample_name}.sam"
done

# Step 2: Convert SAM to BAM, sort, and index
echo "Step 2: Processing SAM files (convert to BAM, sort, and index)..."
for samfile in "${SAMPLEDIR}/result_mapcount/sam/"*.sam; do
    sample_name=$(basename "$samfile" .sam)
    echo "Processing $sample_name..."
    samtools view -bS "$samfile" > "${SAMPLEDIR}/result_mapcount/bam/${sample_name}.bam"
    samtools sort "${SAMPLEDIR}/result_mapcount/bam/${sample_name}.bam" -o "${SAMPLEDIR}/result_mapcount/sorted_bam/${sample_name}_sorted.bam"
    samtools index "${SAMPLEDIR}/result_mapcount/sorted_bam/${sample_name}_sorted.bam"
done

# Step 3: Generate counts using featureCounts
echo "Step 3: Generating counts using featureCounts..."
featureCounts -T "$THREADS" -a "$GTF" -o "${SAMPLEDIR}/result_mapcount/counts/featureCounts.txt" "${SAMPLEDIR}/result_mapcount/sorted_bam/"*.bam

echo "Pipeline completed successfully."
