#!/bin/bash

# Script: normalization_pipeline.sh
# Description: This script normalizes count data for multiple barcode files by calculating geometric means,
#              size factors, and normalized counts.
# Author: Lissa Cruz-Saavedra
# Date: 02-04-2025
# Version: 1.0

# Define the range of barcodes
barcodes=$(seq -w 1 96)

# Input directory for count files
input_dir="input_counts"
output_dir="normalized_counts"

# Ensure input and output directories exist
if [[ ! -d "$input_dir" ]]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi
mkdir -p "$output_dir"

# Function: Calculate geometric means for each gene
calculate_geometric_means() {
    local input_files=("$@")
    awk '
    BEGIN {OFS="\t"}
    {
        # For all count columns, calculate log sum
        for (i=7; i<=NF; i++) { if ($i > 0) { log_sum[$1] += log($i); count[$1]++ } }
    }
    END {
        # Compute geometric means only for genes with counts > 0
        for (gene in log_sum) {
            if (count[gene] > 0) {
                gm = exp(log_sum[gene] / count[gene]);
                print gene, gm
            }
        }
    }
    ' "${input_files[@]}" > "${output_dir}/geometric_means.txt"
}

# Main pipeline
input_files=()
for barcode in $barcodes; do
    input_file="${input_dir}/counts_barcode${barcode}.txt"
    if [[ -f "$input_file" ]]; then
        input_files+=("$input_file")
    else
        echo "Warning: File $input_file not found. Skipping."
    fi
done

if [[ ${#input_files[@]} -eq 0 ]]; then
    echo "Error: No input files found. Exiting."
    exit 1
fi

# Calculate geometric means
calculate_geometric_means "${input_files[@]}"

echo "Geometric means calculated and saved to '${output_dir}/geometric_means.txt'."
