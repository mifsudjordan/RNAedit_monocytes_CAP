#!/bin/bash
set -eu

SECONDS=0

# File paths

vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/hosp_admis/vcf_"
destination_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/all_atoi_strict"


# Initialize an array to store the SRR numbers
declare -a srr_numbers


# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/srr_numbers_hosp.txt"

for srr in "${srr_numbers[@]}"; do

    work_dir=""$vcf_path""$srr""

    cd "$work_dir"

    # transfer to common directory
    cp "$srr"_AtoI_strict.vcf $destination_path

done


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
