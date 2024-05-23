#!/bin/bash
set -eu

SECONDS=0

# Initialize an array to store the SRR numbers
declare -a srr_numbers

# Read the SRR file line by line and saves each number in the array
while IFS= read -r line; do
    srr_numbers+=("$line")
done < "srr_files_one.txt"


for srr in "${srr_numbers[@]}"; do
    echo "Downloading $srr"
    fasterq-dump $srr
    echo "Downloaded $srr successfully."
    gsutil cp *.fastq gs://jorsthesisbucket/
    echo "Running fastqc"
    fastqc *.fastq -o ./fastqc_output

    # 1st pass
    STAR --runThreadN 30 \
         --genomeDir ./index \
         --readFilesIn *.fastq \ 
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix alignment_output1p/ \
         --limitBAMsortRAM 60000000000
    
    echo "Filtering splice junction file"
    cd alignment_output1p
    cat SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > SJ_out_filtered.tab
    cd ..
    
    # 2nd pass
    STAR --runThreadN 30 \
         --sjdbFileChrStartEnd ./alignment_output_1001/SJ_out_filtered.tab \
         --genomeDir ./index \
         --readFilesIn 1001.fastq \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix alignment_output_1001_2p/ \
         --limitBAMsortRAM 60000000000    
    
    #filtering out unaligned reads
    cd alignment_output_1001_2p
    samtools view -F 4 -b Aligned.sortedByCoord.out.bam > aligned.reads.1001.bam

    #indexing bam file
    samtools index aligned.reads.1001.bam

    # Markduplicates

    # SplitNcigar
    gatk SplitNCigarReads -R /home/mifsud_jordan/rawdata/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -I Aligned.sortedByCoord.out.bam -O cigar_output.bam
    
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."