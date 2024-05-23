#!/bin/bash

# see if possible to delete extra file from from bucket before getting the folder to use.
# change name of vcf directory before cp to bucket.

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
done < "srr_numbers.txt"

for srr in "${srr_numbers[@]}"; do
   # Saving directory name of alignment 2nd pass as $output_dir2P
    output_dir2P="aligned_2P${srr}/"

   # Saving name of bam file as variable $aligned_reads
    aligned_reads="aligned.reads.${srr}.bam"
    echo "Hello"
    
   # Get folder from bucket back to instance # 
    gsutil -m cp -r gs://jorsthesisbucket/split_bams/"$output_dir2P" /home/mifsud_jordan/rawdata/

   # ------CD to alignment 2nd pass output directory------ #
    cd "$output_dir2P"

   # Saving output and input file names
    haplo_output="haplo_out_${srr}.vcf"
    recal_rg_output="recal_rg_split_aligned_${srr}.bam"

   # File to which assembled haplotypes will be written
    bamout="bamout_${srr}.bam"
   
   # Running Haplotypecaller
    echo "Running Haplotypecaller"
    $gatk --java-options "-Xmx28g -XX:ParallelGCThreads=8" HaplotypeCaller  \
        -R "$ref_gen" \
        -I "${recal_rg_output}" \
        -O "${haplo_output}" \
        -D "${db_var_snp}" \
        -stand-call-conf 20 \
        -bamout "${bamout}"

    echo "VCF of ${srr} and bamout outputted."

    rm "${recal_rg_output}"
    cd ..
   # ------CD back out to ~./raw_data------ #
   
   # Renaming directory
    vcf_directory="vcf_${srr}"
   
    mv "$output_dir2P" "${vcf_directory}"

   # Copying aligned_bam to bucket
    echo "Transfering vcf directory to bucket"
    gsutil cp -r "$vcf_directory" gs://jorsthesisbucket/vcfs/hosp_admis
   
   # Deleting output2P directory
    rm -r "$vcf_directory"
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
