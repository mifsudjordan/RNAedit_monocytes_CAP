#!/bin/bash

# Program paths

download_location="/mnt/c/Users/jorsm/Downloads/bam_files"

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "srr_numbers.txt"


for srr in "${srr_numbers[@]}"; do
    echo "Indexing $srr"
    cd $download_location/aligned_2P$srr
    samtools index aligned.reads.$srr.bam
done
