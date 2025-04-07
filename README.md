# MGF-seq
Code for analysis of multigene family sequence from Trypanosoma cruzi using nanopore data

# **MGFseq Complete Pipeline: Mapping, Counting, Normalization, and Differential Expression Analysis**

This repository contains a comprehensive pipeline for processing MGF-seq (_Trypanosoma cruzi_ Multigene Family RNA Sequencing) data. The pipeline includes steps for read mapping, BAM file processing, feature counting, normalization, and differential expression analysis. It is designed to handle high-throughput sequencing data efficiently and produce normalized counts and differential expression results for downstream analysis.

# **Features:**

* Read Mapping:
Uses minimap2 to map sequencing reads (FASTQ files) to a reference genome.
Outputs SAM files for further processing.

* BAM File Processing:
Converts SAM files to BAM format using samtools.
Sorts and indexes BAM files for efficient downstream analysis.

* Feature Counting:
Uses featureCounts to count reads mapped to genomic features (e.g., genes) based on an annotation file (GTF format).
Outputs raw counts for each feature.

* Normalization:
Calculates geometric means for each gene across multiple samples.
Normalizes counts using size factors to account for differences in sequencing depth.
The pipeline follows a median of ratios normalization.

* Differential Expression Analysis:
Performs pairwise comparisons between defined groups.
Calculates log fold change (logFC) and p-values for each gene.
Outputs a comprehensive table of differential expression results.

* Error Handling:
Skips missing files with warnings and exits gracefully if no valid input files are found.

# **How It Works:**

**Input:**
FASTQ Files: Raw sequencing reads for each sample.
Reference Genome: A FASTA file for mapping reads.
Annotation File: A GTF file specifying genomic features for counting.
Count Files: Raw count files for each barcode/sample.

**Output:**
SAM Files: Intermediate files from the mapping step.
BAM Files: Processed and sorted BAM files.
Feature Counts: Raw counts for each gene in all samples.
Geometric Means: A file containing geometric means for each gene.
Differential Expression Results: A table of logFC and p-values for all pairwise comparisons.

# **Pipeline Steps:**

* _Mapping Reads:_

Maps reads from FASTQ files to the reference genome using minimap2.
Outputs SAM files to the result_mapcount/sam/ directory.

* _Processing BAM Files:_

Converts SAM files to BAM format using samtools view.
Sorts BAM files using samtools sort.
Indexes BAM files using samtools index.
Outputs sorted BAM files to the result_mapcount/sorted_bam/ directory.

* _Feature Counting:_

Counts reads mapped to genomic features using featureCounts.
Outputs raw counts to result_mapcount/counts/featureCounts.txt.

* _Normalization:_

Calculates geometric means for each gene across all samples.
Normalizes counts using size factors to account for sequencing depth differences.
Outputs normalized counts to normalized_counts/.

* _Differential Expression Analysis:_

Uses an embedded R script to calculate log fold change (logFC) and p-values for each gene between defined groups.
Outputs differential expression results to differential_expression/differential_expression_results.txt.

# Contributing
If you'd like to contribute, please fork the repository and use a feature branch. Pull requests are warmly welcome.

# Licensing
The code in this project is licensed under a CC BY-NC 4.0 license. This means you can use it with proper attribution to original authors, indicating if any modifications were done, only for non commercial applications. (https://creativecommons.org/licenses/by-nc/4.0/)

# Links
If you want to see more of our work you can check out our website: https://www.cestarilab.com/

