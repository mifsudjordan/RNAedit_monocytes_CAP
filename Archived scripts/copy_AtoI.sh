#!/bin/bash
set -eu

SECONDS=0

# File paths

vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/recovery/vcf_"
intersect_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/recovery/intersect_atoi"


# Initialize an array to store the SRR numbers
declare -a srr_numbers


# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/srr_numbers_rec.txt"

for srr in "${srr_numbers[@]}"; do

    work_dir=""$vcf_path""$srr""

    cd "$work_dir"

    # compress
    bgzip -c "$srr"_AtoI.vcf > "$srr"_AtoI.vcf.gz

    # index
    bcftools index "$srr"_AtoI.vcf.gz

    # transfer to common directory
    cp "$srr"_AtoI.vcf $intersect_path
    cp "$srr"_AtoI.vcf.gz $intersect_path
    cp "$srr"_AtoI.vcf.gz.csi $intersect_path

done


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
