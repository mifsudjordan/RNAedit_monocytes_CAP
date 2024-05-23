#!/bin/bash

# Program paths

download_location="/mnt/c/Users/jorsm/Downloads/"

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "srr_numbers.txt"


for srr in "${srr_numbers[@]}"; do
    echo "Downloading $srr"
    mkdir aligned_2P$srr
    gsutil -m cp -r gs://jorsthesisbucket/aligned_2P$srr/* /mnt/c/Users/jorsm/Downloads/aligned_2P$srr
done
