#!/bin/bash
set -eu

SECONDS=0

# File paths

vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/controls/vcf_"
copy_to_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/controls/intersect/"

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "srr_numbers_controls.txt"

for srr in "${srr_numbers[@]}"; do
    
    # copying
    file_to_copy=""$vcf_path""$srr"/"$srr"_CtoU2.vcf"
    cp $file_to_copy $copy_to_path

    # compress
    bgzip -c "$copy_to_path""$srr"_CtoU2.vcf > "$copy_to_path""$srr"_CtoU2.vcf.gz

    # index
    bcftools index "$copy_to_path""$srr"_CtoU2.vcf.gz

done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

