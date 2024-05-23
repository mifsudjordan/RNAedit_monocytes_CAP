# RNAedit_monocytes_CAP
Contains scripts used for my bioinformatics masters thesis entitled "Investigating RNA editing in human primary monocytes during community-acquired pneumonia".
This repository contains a set of scripts used for discovering RNA editing events in RNA-seq samples. The scripts are organized to form a pipeline for processing and analyzing RNA-seq data.

## Prerequisites
Unix-based operating system (Linux/MacOS)
Bash shell
Python 3.x
R
Required Python and R libraries (listed in individual scripts)

## Pipeline Overview
The pipeline consists of the following main steps:

1. Download and alignment of RNA-seq data.
2. Splitting of reads and base recalibration.
3. Variant calling.
4. Variant filtering and counting.
5. Compression, copying, intersection, and annotation of VCF files.
6. Processing of intersection results.
7. Statistical analysis of results.

## Main Scripts
### 1. 1_srr_download_align.sh
This script performs several tasks:
- Download RNA-seq data from the SRA (Sequence Read Archive).
- Run quality control on the downloaded data using FastQC.
- Align reads to a reference genome using STAR.
- Sort and filter the aligned reads using Samtools.
- Upload results to a Google Cloud Storage bucket.

What it works: 
The script reads SRR numbers from srr_numbers.txt and stores them in an array.
A processing loop performs the following for every SRR number:
- Download the FASTQ file using fasterq-dump.
- Run FastQC on the downloaded FASTQ file.
- First pass alignment with STAR to generate an unsorted BAM file and splice junctions.
- Filter splice junctions based on specified criteria.
- Second pass alignment with STAR using filtered splice junctions.
- Sort the BAM file using Samtools.
- Filter out unaligned reads to generate the final BAM file.
- Upload results to the specified Google Cloud Storage bucket.
- Clean up intermediate files to save space.

Note:
Ensure you have the required tools installed and accessible in your system's PATH.
Prepare a file named srr_numbers.txt containing the SRR numbers to process, one per line.

Example srr_numbers.txt:

SRR123456

SRR123457

SRR123458
etc.

Also ensure the paths to gatk, ref_gen, and db_var_snp are correctly set in the script.
Modify the Google Cloud Storage bucket paths (gs://jorsthesisbucket/) as needed.


Dependencies:

- GATK
- Fasterq-dump
- FastQC
- STAR
- Samtools
- gsutil (for Google Cloud Storage operations)

### 2. 2_splitreads_baseRecalibration.sh
This script performs several tasks to prepare aligned reads for variant calling:
- Splits reads that span introns in RNA-seq data.
- Adds read group information to the BAM file.
- Performs base quality score recalibration using GATK.
- Uploads the processed files to a Google Cloud Storage bucket.

How it works:
- Reads SRR numbers from srr_numbers_controls3.txt and stores them in an array.

For each SRR number in the array:
- Download alignment results from a Google Cloud Storage bucket.
- Change to the alignment output directory.
- Define file names for input and output BAM files.
- Run SplitNCigarReads: Splits reads that span introns to make the data compatible with GATK tools.
- Add Read Groups: Adds metadata to the BAM file using AddOrReplaceReadGroups.
- Base Recalibration: Runs BaseRecalibrator to create a recalibration table.
- Applies the recalibration using ApplyBQSR.
- Clean Up: Removes intermediate files to save space.
- Upload Results: Transfers the processed BAM file back to the Google Cloud Storage bucket.
- Clean Up: Removes the local copy of the alignment directory.
- Completion: Prints the total execution time.

*to be continued*

## Custom Utilities
The custom_utilities directory contains various helper scripts used by the main scripts or for small alterations to the output files. These include scripts for downloading data, finding specific variants, and generating statistical tables.

## Archived Scripts
The archived scripts directory contains previous versions of the main scripts and drafts. These scripts are preserved for reference and may contain useful code snippets or alternative approaches to the tasks performed by the main scripts.

Usage
To run the pipeline, execute the main scripts in the order listed above and change paths and file names according to your needs. Ensure that all dependencies and prerequisites are installed and configured correctly.
