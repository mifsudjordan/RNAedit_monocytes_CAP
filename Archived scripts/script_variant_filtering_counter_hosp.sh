#!/bin/bash
set -eu

SECONDS=0

# File paths

ref="/mnt/c/Users/jorsm/OneDrive/masters_thesis/references_indexes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
vcf_path="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/hosp_admis/vcf_"
output="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/hosp_admis/variant_stats_hosp3.txt"
gatk="/root/gatk-4.4.0.0/gatk"
check_REDI="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/check_REDI.py"

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

    work_dir=""$vcf_path""$srr""

    cd "$work_dir"

    echo "$srr" >> $output
    echo "No. unfiltered variants inc. INDELS" >> $output
    grep AC= haplo_out_"$srr".vcf | wc -l >> $output
    
    # Filtering out indels and variants in coordinates of common SNPs 
    # according to dbSNP database.

    $gatk SelectVariants \
        -V haplo_out_"$srr".vcf \
        -ids . \
        -select-type SNP \
        -O "$srr"_uncommon_snps.vcf

    
    echo "No. of uncommon SNPs" >> $output
    grep AC= "$srr"_uncommon_snps.vcf | wc -l >> $output

    # QC filtration:
    $gatk VariantFiltration \
        -R "$ref" \
        -V "$srr"_uncommon_snps.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O "$srr"_uncommonsnps_qc.vcf

    # select only the passes...

    $gatk SelectVariants \
        -V "$srr"_uncommonsnps_qc.vcf \
        --exclude-filtered \
        -O "$srr"_uncommon_snps_passes.vcf

    echo "No. of uncommon SNPs that passed QC" >> $output
    grep AC= "$srr"_uncommon_snps_passes.vcf | wc -l >> $output

    # Selecting specific variants:

    echo "Selecting variants matching RNA editing types..."

    # A to I || C to U
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} ($4 == "G" && $5 == "T") || ($4 == "C" && $5 == "A") || ($4 == "A" && $5 == "G") || ($4 == "T" && $5 == "C")' "$srr"_uncommon_snps_passes.vcf > "$srr"_allRNAedits.vcf

    echo "No. of A to I || C to U" >> $output
    grep AC= "$srr"_allRNAedits.vcf | wc -l >> $output

    # A to I
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} ($4 == "T" && $5 == "C") || ($4 == "A" && $5 == "G")' "$srr"_uncommon_snps_passes.vcf > "$srr"_AtoI.vcf

    echo "No. of A to I" >> $output
    grep AC= "$srr"_AtoI.vcf | wc -l >> $output

    # C to U
    awk 'BEGIN {OFS="\t"} /^#/ {print; next} ($1 ~ /^[0-9XYM]+$/) && (($4 == "G" && $5 == "T") || ($4 == "C" && $5 == "A"))' "$srr"_uncommon_snps_passes.vcf > "$srr"_CtoU2.vcf
    echo "No. of C to U on chromosomes only:" >> $output
    grep AC= "$srr"_CtoU2.vcf | wc -l >> $output

    echo "Counting matching variants to REDI portal database..."

    $check_REDI "$srr"_AtoI.vcf >> $output

done
    echo "-----------------------" >> $output

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

    # Annotation by gene
    # annovar/table_annovar.pl "$srr"_allRNAedits.vcf annovar/humandb/ -buildver hg38 -out "$srr"_allRNAedits_ann --thread 4 -remove -protocol refGene -operation g -nastring . -vcfinput -polish
