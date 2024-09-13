# RNAedit_monocytes_CAP
Contains scripts used for my bioinformatics masters thesis entitled "Investigating RNA editing in human primary monocytes during community-acquired pneumonia".
This repository contains a set of scripts used for discovering RNA editing events in RNA-seq samples. The scripts are organized to form a pipeline for processing and analyzing RNA-seq data.

## Prerequisites
- Unix-based operating system (Linux/MacOS)
- Bash shell
- Python 3.x
- R
- Required Python and R libraries (listed in individual scripts)
- R studio

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

How it works: 
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

- `SRAtoolkit fasterq-dump`
- `FastQC`
- `STAR`
- `Samtools`
- `gsutil` (for Google Cloud Storage operations)

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

Dependencies:
- `gatk`
- `gsutil` (for Google Cloud Storage operations)

### 3. 3_variant_calling.sh
This script performs one main task:
- Performs variant calling on the output of the previous script

How it works:
- Reads SRR numbers from srr_numbers_controls3.txt and stores them in an array.

For each SRR number in the array:
- Downloads folder containing output of previous script from a Google Cloud Storage bucket.
- Runs HaplotypeCaller to output vcf files from the processed BAM files.
- Clean Up: Removes intermediate files to save space.
- Upload vcf directory back to the Google Cloud Storage bucket.
- Clean Up: Removes the local copy of the alignment directory.
- Completion: Prints the total execution time.

Dependencies:
- `gatk`
- `gsutil` (for Google Cloud Storage operations)

The following scripts were run locally after downloading all vcf output files from cloud storage.

### 4. 4_variant_filtering_counter.sh

This script performs several tasks:
- Variant quality control
- Variant filtration
- Counts variants showing RNA editing patterns
- Outputs new vcfs with filtered variants
- Counts potential RNA editing variants
- Outputs RNA editing counts results for every sample

How it works:
- Reads SRR numbers from srr_numbers_controls3.txt and stores them in an array.

For each SRR number in the array:
- The working directory name was amended as per sample file
- gatk VariantFiltration marks each variant as PASS or not
- gatk SelectVariants outputs vcf files with only PASSED variants
- Counts total potential RNA edits
- Counts potential A to I RNA edits
- Counts potential C to U RNA edits
- Counts A to I edits that match the REDI portal database v2.
- Outputs all counts in a results file.

Dependencies:
- `gatk`
- `check_REDI.py` (from the Custom Utilities directory)

Note: The results file outputted by the script contains the counts of RNA edits for every sample.
The file can be converted to a tab-delimited txt or csv file that can easily interpreted by an ETL tool such as Power Query.

### 5. 5_vcf_compressor_copier_intersect_annotator.sh

This script performs:
- vcf file compression
- vcf file indexing
- Transfers files from a list of directories to a common directory
- vcf file intersection
- vcf file format conversion
- variant annotation

How it works:
- Reads SRR numbers from srr_numbers_controls3.txt and stores them in an array.

For each SRR number in the array:
- Visits the directory belonging to the sample
- Bgzips the vcf file
- Indexes the compressed vcf file using bcftools index
- Copies files to a common directory

  Once ready the script transfers to the common directory containing all compressed indexed vcfs and
- Intersects the vcf files using bcftools isec
- Converts the vcf file format to one required by annotation program annovar, using a custom python program intersect_to_annovar.py
- Annotates (gene annotation) the variants in this modified vcf using annovar  

Dependencies:
- `bgzip`
- `bcftools`
- `annovar`
- `intersect_to_annovar.py` (from the Custom Utilities directory)

### 6. 6_intersect_processor2.py

This script performs:
- Processing of intersect file outputted by bcftools isec
- Outputs a tab-delimited csv file with all annotations and statistics per RNA edit

How it works:
The program uses a loads a text file containing the group / time-point to which every sample belongs to (acute/recovery/control) sorted in the same order as the vcf files have been processed by bcftools isec command. Then it loads the intersect file (bcftools isec output) which contains a string of 1s and 0s per RNA edit indicating the presence or absence of that RNA edit in the samples. The order of 1s and 0s is the same as the vcf file list processed by bcftools (e.g. sorted in ascending order by filename. Hence if the string is 00010010, then the RNA edit is found in the 4th and 7th vcf file analysed by bcftools. 

For each string of 1s and 0s per RNA edit:
- Identifies the group (acute/recovery/control) in which that edit was detected
- Counts the occurances and finds the percentage occurance of that edit in the whole study
- Counts the occurances and finds the percentage occurance of that edit per group / time-point.

Once ready the program outputs a tab-delimited csv file with all annotated RNA edits and statistics 

Dependencies:
- a .txt file with the time-points/group of each sample
- intersect file output

### 7. 7_thesis_rscript_strict.R / 7_thesis_rscript_norm.R

These scripts performs similar data exploration and analysis on the outputs of the previous scripts and metadata of each sample. 

Datasets with counts of RNA edits per chromosome per sample were combined with sample metadata and four cytokine-release values of monocytes upon LPS-stimulation (TNF-alpha, IL-1b, Il-6 and Il-10) in preparation for a multi-omics analysis. Statistical analysis of the RNA edit counts and cytokine-release data was carried out using these two custom R scripts for the count data derived using the normal and strict variant filtering approaches. Quantitative analysis included defining descriptive statistics, correlation analysis, principal component analysis (PCA) and hypothesis testing. After considering the types of data available and their measures of dispersion, for most hypothesis tests, the Negative Binomial Test was used. This is a generalised linear model test, often used for count data that shows over-dispersion. 

Alignment data was used to count reads that mapped on different genes using Subread’s featureCounts on the command-line (not in R script). Genes with zero expression were removed from all samples. Gene names were obtained using the org.Hs.eg.db annotation R package from Bioconductor. A design matrix was constructed to incorporate experimental conditions, specifying the timepoint each sample belonged to. This allowed for the modelling of differential expression across groups.
Differential expression analysis was performed using Bioconductor’s DESeq2 package. This included estimating size factors to account for sequencing depth differences and fitting a negative binomial model to assess differential expression. Significance testing was conducted to identify genes with expression changes between experimental conditions, with adjustments for multiple testing. PCA and dispersion estimates were used to visualize the results and assess the quality of the differential expression analysis.
Code in these R-scripts generate bar plots to illustrate significant log2fold changes between experimental conditions for specific genes of interest, focusing on RNA editing genes. Normalized gene expression counts were extracted from the differential expression analysis and were merged with RNA editing count data and cytokine-release data for a multiomics analysis. Correlation analysis was conducted using normalised RNA-editing gene counts, RNA-edit counts per chromosome and cytokine-release data, while distinguishing between different timepoints. PCA was conducted on cytokine-release data and normalised counts of RNA-editing genes. This comprehensive approach enabled the identification of associations between RNA editing counts, cytokine-release capacities and gene expression, particularly focusing on RNA editing genes, across different experimental conditions.

Code in these R scripts were used to process the output of the intersection file to elucidate differential occurrence of specific edits between the different time-points. Fisher’s exact test was used to check the statistical significance of such differences. P-values were adjusted for multiple hypothesis testing using the Benjamin-Hochberg correction. Mann-Witney U tests were carried out to find if cytokine-release values of such samples were significantly different from those of other samples within the same group that did not have the same RNA edit. 

Dependencies:
- `tidyverse`
- `psych`
- `readr`
- `ggplot2`
- `tidyr`
- `DESeq2`
- `apeglm`
- `biomaRt`
- `gplots`
- `gtools`
- `ggrepel`
- `Hmisc`
- `MASS`

## Custom Utilities
The custom_utilities directory contains various helper scripts used by the main scripts or for small alterations to the output files. These include scripts for simple ETL tasks, finding sample names containing RNA edits, comparing vcf records to databases such as REDIportal database and generating descriptive statistical tables for each RNA edit.

## Archived Scripts
The archived scripts directory contains previous versions of the main scripts and drafts. These scripts are preserved for reference and may contain useful code snippets or alternative approaches to the tasks performed by the main scripts.

# Usage
To run the pipeline, execute the main scripts in the order listed above and change paths and file names according to your needs. Ensure that all dependencies and prerequisites are installed and configured correctly.
For the R script, load the scripts in R studio, manaually executing the code blocks of the analyses you're interested to perform. It is not suggested that the R scripts are executed as a whole all at once. Many parts of the scripts have been used for data exploratory purposes.   
