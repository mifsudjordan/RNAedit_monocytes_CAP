#!/bin/bash
set -eu

SECONDS=0

# File paths

ref="/mnt/c/Users/jorsm/OneDrive/masters_thesis/references_indexes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/hosp_admis/vcf_"
output="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/hosp_admis/variant_stats_acute_norm.txt"
gatk="/root/gatk-4.4.0.0/gatk"
check_REDI="/mnt/c/Users/jorsm/OneDrive/masters_thesis/scripts/from_onedrive_folder/vcf_scripts/check_REDI.py"

# Initialize an array to store the SRR numbers
declare -a srr_numbers

touch $output
echo $(date) >> $output

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    if [ -n "$line" ]; then
        srr_numbers+=("$line")
    fi
done < "srr_numbers_hosp.txt"

for srr in "${srr_numbers[@]}"; do
    echo $srr >> $output
    
    work_dir=""$vcf_path""$srr""

    cd "$work_dir"

    # QC normal filtering
    $gatk VariantFiltration \
        -R "$ref" \
        -V "$srr"_uncommon_snps.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O "$srr"_uncommonsnps_qc_norm.vcf

    # QC strict filtering:
#    $gatk VariantFiltration \
#        -R "$ref" \
#        -V "$srr"_uncommon_snps.vcf \
#        -filter "QD < 2.0" --filter-name "QD2" \
#        -filter "DP < 10" --filter-name "DP10*" \
#        -filter "HRun > 5" --filter-name "HRun5*" \
#        -filter "MQ < 40.0" --filter-name "MQ40" \
#        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#        -filter "AD[1] < 3" --filter-name "AltAlleleSupportLow*" \
#        -filter "AD[1]/(AD[0]+AD[1]) > 0.9" --filter-name "AltAlleleSupportHigh*" \
#        -O "$srr"_uncommonsnps_qc_new.vcf

    # select only the passes...

    $gatk SelectVariants \
        -V "$srr"_uncommonsnps_qc_norm.vcf \
        --exclude-filtered \
        -O "$srr"_uncommon_snps_passes_norm.vcf

    echo "No. of uncommon SNPs that passed QC" >> $output
    grep AC= "$srr"_uncommon_snps_passes_norm.vcf | wc -l >> $output

    # Selecting specific variants:

    echo "Selecting variants matching RNA editing types..."

    # A to I || C to U
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} ($4 == "G" && $5 == "T") || ($4 == "C" && $5 == "A") || ($4 == "A" && $5 == "G") || ($4 == "T" && $5 == "C")' "$srr"_uncommon_snps_passes_norm.vcf > "$srr"_allRNAedits_norm.vcf

    echo "No. of A to I || C to U" >> $output
    grep AC= "$srr"_allRNAedits_norm.vcf | wc -l >> $output

    # A to I
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} ($4 == "T" && $5 == "C") || ($4 == "A" && $5 == "G")' "$srr"_uncommon_snps_passes_norm.vcf > "$srr"_AtoI_norm.vcf

    echo "No. of A to I:" >> $output
    grep AC= "$srr"_AtoI_norm.vcf | wc -l >> $output

    # C to U
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} (($1 ~ /^[0-9XY]+$/) || ($1 == "MT")) && (($4 == "G" && $5 == "T") || ($4 == "C" && $5 == "A"))' "$srr"_uncommon_snps_passes_norm.vcf > "$srr"_CtoUnorm.vcf
    echo "No. of C to U:" >> $output
    grep AC= "$srr"_CtoUnorm.vcf | wc -l >> $output

    echo "Counting matching variants to REDI portal database..."

    $check_REDI "$srr"_AtoI_norm.vcf >> $output

done
    echo "-----------------------" >> $output

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

    # Annotation by gene
    # annovar/table_annovar.pl "$srr"_allRNAedits.vcf annovar/humandb/ -buildver hg38 -out "$srr"_allRNAedits_ann --thread 4 -remove -protocol refGene -operation g -nastring . -vcfinput -polish
