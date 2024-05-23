#!/bin/bash
set -eu

SECONDS=0

# File paths

vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/recovery/vcf_"
intersect_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/recovery/intersect_atoi_recovery_norm"
srr_numbers_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/"
intersect_to_annovar="/mnt/c/Users/jorsm/OneDrive/masters_thesis/scripts/custom_utilities/intersect_to_annovar.py"


# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "$srr_numbers_path""srr_numbers_rec.txt"

for srr in "${srr_numbers[@]}"; do

    work_dir=""$vcf_path""$srr""

    cd "$work_dir"

    # compress
    bgzip -c "$srr"_AtoI_norm.vcf > "$srr"_AtoI_norm.vcf.gz

    # index
    bcftools index "$srr"_AtoI_norm.vcf.gz

    # transfer to common directory
    cp "$srr"_AtoI_norm.vcf $intersect_path
    cp "$srr"_AtoI_norm.vcf.gz $intersect_path
    cp "$srr"_AtoI_norm.vcf.gz.csi $intersect_path

done

cd "$intersect_path"

# intersect all .vcf.gz files in the pwd
bcftools isec -n +2 *.vcf.gz > intersect_atoi_all.vcf

#converting intersect file to format required by annovar
$intersect_to_annovar intersect_atoi_all.vcf

# annotation of processed intersect file
annovar/table_annovar.pl intersect_annovar_ctou_norm.vcf annovar/humandb/ \
                         -buildver hg38 \
                         -out annovar_output/intersect_ctou_ann \
                         -remove \
                         -protocol refGene \
                         -operation g \
                         --otherinfo


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
