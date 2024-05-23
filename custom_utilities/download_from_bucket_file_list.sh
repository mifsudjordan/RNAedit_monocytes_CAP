#!/bin/bash

# Local directory where you want to download the files
local_directory="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/hosp_admis"

# Bucket name
bucket_name="gs://jorsthesisbucket/vcfs/hosp_admis"

# Loop through the list of filenames in left.txt and download them
while IFS= read -r filename
do
  gsutil -m cp -r "$bucket_name"/"$filename" "$local_directory"
done < left.txt