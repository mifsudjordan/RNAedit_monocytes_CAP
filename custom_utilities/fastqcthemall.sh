#!/bin/bash

output_dir="/mnt/c/Users/jorsm/OneDrive/masters_thesis/rawdata/fastqc_output"
# creates output directory if not present
mkdir -p "$output_dir"

echo "Decompressing contents of directory and outputting in"
echo $output_dir


file_count=0
for fastq_file in *.fastq; do
    ((file_count = file_count + 1))
done
echo "There are $file_count fastq files in directory"

counter=0
for fastq_file in *.fastq; do
    # "file" will store the filenames without the .fastq extension 
    fastqc $fastq_file -o $output_dir
    ((counter = counter + 1))
    echo "$counter out of $file_count processed."
    echo "Processed $fastq_file"
done

echo "Processed all files."