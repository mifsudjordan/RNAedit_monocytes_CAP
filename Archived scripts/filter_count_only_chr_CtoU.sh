#!/bin/bash
set -eu

SECONDS=0

# File paths

vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/controls/vcf_"
output="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/controls/count_CtoUchr_only2.txt"


# Initialize an array to store the SRR numbers
declare -a srr_numbers

touch $output
echo $(date) >> $output

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "srr_numbers_controls.txt"

for srr in "${srr_numbers[@]}"; do

    work_dir=""$vcf_path""$srr""

    cd "$work_dir"

    echo "$srr" >> $output
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} ($1 ~ /^[0-9XYM]+$/) && (($4 == "G" && $5 == "T") || ($4 == "C" && $5 == "A"))' "$srr"_uncommon_snps_passes.vcf > "$srr"_CtoU2.vcf
    echo "No. of C to U on chromosomes only:" >> $output
    grep AC= "$srr"_CtoU2.vcf | wc -l >> $output

done
    echo "-----------------------" >> $output

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

    # Annotation by gene
    # annovar/table_annovar.pl "$srr"_allRNAedits.vcf annovar/humandb/ -buildver hg38 -out "$srr"_allRNAedits_ann --thread 4 -remove -protocol refGene -operation g -nastring . -vcfinput -polish
