#!/bin/bash
set -eu

SECONDS=0

# File paths
intersect_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/all_ctou_strict_compressed/"
intersect_to_annovar="/mnt/c/Users/jorsm/OneDrive/masters_thesis/scripts/RNAedit_monocytes_CAP/custom_utilities/intersect_to_annovar.py"

# annotation of processed intersect file
annovar/table_annovar.pl intersect_annovar_atoi_strict.vcf annovar/humandb/ \
                         -buildver hg38 \
                         -out annovar_output/intersect_atoi_strict_ann \
                         -remove \
                         -protocol refGene \
                         -operation g \
                         --otherinfo


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
