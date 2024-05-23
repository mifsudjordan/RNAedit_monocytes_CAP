#!/bin/bash

# Program paths
download_location="/mnt/c/Users/jorsm/Downloads/bam_files/"
paste_to="/mnt/c/Users/jorsm/Downloads/just_the_bams/"

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and save each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "/mnt/c/Users/jorsm/Downloads/srr_numbers.txt"

# Check if the destination directory exists, and if not, create it
if [ ! -d "$paste_to" ]; then
    mkdir -p "$paste_to"
fi

for srr in "${srr_numbers[@]}"; do
    echo "Moving $srr"
    cp -r "$download_location"aligned_2P$srr/*.bam "$paste_to"aligned.reads.$srr.bam
    cp -r "$download_location"aligned_2P$srr/aligned.reads.$srr.bam.bai "$paste_to"aligned.reads.$srr.bam.bai
    rm "$download_location"aligned_2P$srr/aligned.reads.$srr.bam
done

echo "All files moved."