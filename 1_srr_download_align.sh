#!/bin/bash
set -eu # Script will halt upon encountering an error
        # including undefined variables. 

SECONDS=0

# ulimit -n 10000 # Sets the limit to open file descriptors
                  # beyond the usual defaults of 1024 or 4096.
                  # I hashed it out as it was probably overkill 
                  # for this use case.

# Program paths
gatk="/home/mifsud_jordan/gatk-4.4.0.0/gatk"
ref_gen="/home/mifsud_jordan/rawdata/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
db_var_snp="/home/mifsud_jordan/rawdata/00-common_all.vcf"
programs=("fasterq-dump" "fastqc" "STAR" "samtools")

# Checking if all programs are working:
for program in "${programs[@]}"; do
    if command -v "$program" &>/dev/null; then
        echo "Program '$program' is installed."
        version="$("$program" --version)"
        echo "Version: $version"
        echo "---------------------------------------------"
    else
        echo "*********************************************"
        echo "Program '$program' is not installed or not in the system's PATH."
        echo "*********************************************"
    fi
done
$gatk --version

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array.
while IFS= read -r line; do # Reads the input file line by line preserving
                            # whitespace and backslashes.
    if [ -n "$line" ]; then # prevents blank lines from being added to array
        srr_numbers+=("$line") 
    fi
done < "srr_numbers.txt"


for srr in "${srr_numbers[@]}"; do
   # Downloading fastq file in pwd (~/rawdata)
    echo "Downloading $srr"
    fasterq-dump $srr
    echo "Downloaded $srr successfully."
   
   # Saving name of fastq file as variable $downloaded_file
    downloaded_file="${srr}.fastq"
   
   # Copying fastq file to bucket
    gsutil cp "./$downloaded_file" gs://jorsthesisbucket/

   # Running fastqc on fastq file
    echo "Running fastqc."
    fastqc "$downloaded_file" -o ./fastqc_output
   
   # Saving fastqc name as variable $fastqc_file
    fastqc_file="${srr}_fastqc.html"

   # Copying fastqc output to bucket
    gsutil cp "./fastqc_output/$fastqc_file" gs://jorsthesisbucket/fastqc_output/

   # Saving alignment output directory name in variable $output_dir
    output_dir="aligned_${srr}/"

   # Alignment 1st pass
   # Run STAR to generate an unsorted BAM file
    STAR --runThreadN 20 \
         --genomeDir ./index \
         --readFilesIn "$downloaded_file" \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix "$output_dir" 


   # ------CD to alignment output directory------ #
    cd "$output_dir"

   # Filtering splice junction file
   # Removing Likely false positive splice junctions such as
   # junctions supported by very few reads (e.g. ≤ 2 reads), 
   # non-canonical junctions (The 0 value in column5 from SJ.out.tab denotes non-canonical junctions), 
   # and annotated junctions (The 1 value in column6 from SJ.out.tab denotes annotated junctions). 
   # The annotated junctions should be filtered out as they are already included in the gene annotation file.
    cat SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ_out_filtered.tab
    cd ..
   # ------CD back out to ~/rawdata-------------- #

   # Saving path of filtered splice junction file
    SJ_filtered="./${output_dir}SJ_out_filtered.tab"

   # Saving directory name of alignment 2nd pass as $output_dir2P
    output_dir2P="aligned_2P${srr}/"
    
    # Alignmnet 2nd pass
    STAR --runThreadN 20 \
         --sjdbFileChrStartEnd "$SJ_filtered" \
         --genomeDir ./index \
         --readFilesIn "$downloaded_file" \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix "$output_dir2P" 
    echo "2nd pass alignment completed on $srr"

# Deleting 1st pass alignment output directory
    rm -r "$output_dir"
    rm "$downloaded_file"
    echo "Removed first pass bam and fastq file."

# Sort the unsorted BAM file using Samtools
    echo "Sorting bam file."
    samtools sort -@ 20 -o "${output_dir2P}sorted_output.bam" "${output_dir2P}Aligned.out.bam"

   # ------CD to alignment 2nd pass output directory------ #
    cd "$output_dir2P"
    
   # Saving name of bam file as variable $aligned_reads
    aligned_reads="aligned.reads.${srr}.bam"

   # Removing aligned (unsorted) bam file.
    rm Aligned.out.bam

   # Filtering out unaligned reads
    echo "Filtering out unaligned reads"
    samtools view -F 4 -b sorted_output.bam > "$aligned_reads"
    echo "${aligned_reads} outputted."
  
   # Removing mixed bam file with unaligned reads
    rm sorted_output.bam
   
    cd ..
   # ------CD back out to ~./raw_data------ #

   # Copying aligned_bam to bucket
    echo "Transfering bam to bucket"
    gsutil cp -r "$output_dir2P" gs://jorsthesisbucket/

   # Deleting $output_dir2P directory before restarting loop
    rm -r "$output_dir2P"

done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
