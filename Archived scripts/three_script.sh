#!/bin/bash
set -eu

SECONDS=0

ulimit -n 10000

# Program paths
gatk="/home/mifsud_jordan/gatk-4.4.0.0/gatk"
ref_gen="/home/mifsud_jordan/rawdata/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
db_var_snp="/home/mifsud_jordan/rawdata/00-common_all.vcf"

$gatk --version

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "srr_numbers2.txt"

for srr in "${srr_numbers[@]}"; do
   # Saving directory name of alignment 2nd pass as $output_dir2P
    output_dir2P="aligned_2P${srr}/"
   
   # Get folder from bucket back to instance # 
    gsutil -m cp -r gs://jorsthesisbucket/"$output_dir2P" /home/mifsud_jordan/rawdata/

   # ------CD to alignment 2nd pass output directory------ #
    cd "$output_dir2P"

   # Saving name of bam file as variable $aligned_reads
    aligned_reads="aligned.reads.${srr}.bam"

   # Saving output file name of splitncigarreads as $split_aligned_reads
    split_aligned_reads="split_aligned_${srr}.bam"
   
   # Running SplitNCigarReads
    echo "Running splitncigar"
    $gatk SplitNCigarReads \
        -R "$ref_gen" \
        -I "$aligned_reads" \
        -O "$split_aligned_reads" 
   
   # Removing aligned_reads file
    rm "${aligned_reads}"

   # Saving name of bam file with read groups as variable $rg_output
   # Adding metadata lines in bam file
    rg_output="rg_${split_aligned_reads}"
    $gatk AddOrReplaceReadGroups \
        I="$split_aligned_reads" \
        O="$rg_output" \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=Hosp_admis \
        RGSM="${srr}"

   # Deleting split_aligned_reads bam file
    rm "${split_aligned_reads}"

   # Creating base recalibration data table
    echo "$split_aligned_reads processing by BaseRecalibrator"
    $gatk BaseRecalibrator \
        -I "$rg_output" \
        -R "$ref_gen" \
        --known-sites "$db_var_snp" \
        -O recal_data.table

   # Saving recalibrated bam file name as variable $recal_rg_output
   # Recalibrating base quality scores of bam file.
    recal_rg_output="recal_${rg_output}"
    $gatk ApplyBQSR \
        -R "$ref_gen" \
        -I "$rg_output" \
        --bqsr-recal-file recal_data.table \
        -O "${recal_rg_output}"
        
   # Deleting un-calibrated bam file.
    rm "${rg_output}"
    cd ..
   # ------CD back out to ~./raw_data------ #

   # Copying aligned_bam to bucket
    echo "Transfering bam to bucket"
    gsutil cp -r "$output_dir2P" gs://jorsthesisbucket/split_bams/
   
   # Deleting output2P directory
    rm -r "$output_dir2P"
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
