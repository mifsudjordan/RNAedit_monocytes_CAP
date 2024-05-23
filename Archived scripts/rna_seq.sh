#!/bin/bash

SECONDS=0

echo "Thesis Masters pipeline - Jordan Mifsud"

# going to working directory
cd /mnt/c/Users/jorsm/OneDrive/masters_thesis/rawdata
# transferring compressed files to "compressed_raw" directory
mkdir compressed_raw
mv * compressed_raw/
cd /compressed_raw

# Running decompress.sh to decompress contents of pwd 
# and put them in parent directory
chmod +x decompress.sh
./decompress.sh

# go up one folder
cd ..

# Running fastqc on all files in /rawdata
chmod +x fastqcthemall.sh
./fastqcthemall

###I did not do the following###
# Trimming
java -jar /mnt/wslg/distro/usr/share/java/trimmomatic-0.39.jar SE ENCFF493KQW.fastq ENCFF493KQW_trimmed.fastq ILLUMINACLIP:/mnt/wslg/distro/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10

echo "Trimming completed."

# Running fastqc again on trimmed fastq file
fastqc ENCFF493KQW_trimmed.fastq -o fastqc_output/

echo "Unzipping human reference genome and creating a fasta file containing only"
echo "chromosomes 21 and 22"
###I did not do the above###

# Unzipping file
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
# A soft-masked reference genome refers to a genome sequence in which repetitive and low-complexity 
# regions are marked or "masked" to reduce their influence during sequence analysis. This masking 
# involves replacing or converting repetitive and low-complexity sequences into a uniform representation, 
# often using lowercase letters or other symbols.
# The masking of repetitive and low-complexity regions is done to mitigate potential biases and artifacts 
# that these regions can introduce in various sequence analysis tasks, such as genome alignment, variant calling,
# and gene expression quantification.

echo "Indexing genome."
# Indexing genome
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
# The above command generates an index file for the reference genome in FASTA format, 
# which allows for efficient retrieval of sequences from the reference during various 
# analysis steps. An index file with the .fai extension (e.g., GRCh38_latest_genomic.fna.fai) 
# will be created in the same directory as the reference genome.

# create the dict file that will be used by gatk
# go to the directory were the reference genome is
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa


# decompressing annotation file
gunzip GRCh38_latest_genomic.gff.gz

# Indexing the reference genome file
STAR --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.110.gtf --genomeSAindexNbases 12 --sjdbOverhang 49 --runThreadN 14
# if this does not work:
# try to generate indexes for half of the genome.

echo "Aligning"

# To run alignment of all fastq files in a directory, use the bash script alignment_script.sh
# 1st pass
STAR --runThreadN 30 --genomeDir ./index --readFilesIn 2007.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignment_output/ --limitBAMsortRAM 60000000000

#2nd pass
STAR --runThreadN 30 --sjdbFileChrStartEnd ./alignment_output_1001/SJ_out_filtered.tab --genomeDir ./index --readFilesIn 1001.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignment_output_1001_2p/ --limitBAMsortRAM 60000000000

# checked some stats on each of the aligned sorted files.
# See flagstat output in flagstat_statistics.txt
samtools flagstat Aligned.sortedByCoord.out.bam

# There were no unaligned reads so the following step was skipped.
# If required, discard unaligned reads using samtools
samtools view -F 4 -b aligned_readsAligned.sortedByCoord.out.bam > aligned_reads_only.bam

# Creating bam index files (.bam.bai) using samtools index
samtools index aligned.reads.1001.bam

# Running splitncigarreads
gatk SplitNCigarReads -R /mnt/c/Users/jorsm/OneDrive/masters_thesis/references_indexes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -I Aligned.sortedByCoord.out.bam -O cigar_outpu
t.bam

# or this:
gatk SplitNCigarReads -R /home/mifsud_jordan/rawdata/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -I Aligned.sortedByCoord.out.bam -O cigar_output.bam


# sorting bed and bam files:
sort -k1,1V -k2,2n -k3,3n -o chr21_22_annotation_sorted.bed chr21_22_annotation.bed
samtools sort -o sorted_aligned_reads_only.bam aligned_reads_only.bam

# modified names of chromosomes on the bam file to match annotation file
samtools view -h sorted_aligned_reads_only.bam | sed 's/NC_000021\.9/chr21/g; s/NC_000022\.11/chr22/g' | samtools view -b -o modified_aligned_reads.bam

# sorting, compressing and indexing gtf annotation file

sort -k1,1 -k4,4n -o Homo_sapiens.GRCh38.110.sorted.gtf Homo_sapiens.GRCh38.110.gtf
bgzip -k Homo_sapiens.GRCh38.110.sorted.gtf
tabix -p gff Homo_sapiens.GRCh38.110.sorted.gtf.gz

#trying to use subread's featureCounts
./featureCounts -a /mnt/c/Users/jorsm/OneDrive/masters_thesis/references_indexes/Homo_sapiens.GRCh38.110.sorted.gtf.gz \
  -F GTF -g gene_id -o count.out -T 14 \
  /mnt/c/Users/jorsm/Downloads/just_the_bams/*.bam

echo "Estimating gene coverage for each gene."
bedtools coverage -a chr21_22_annotation_sorted.bed -b modified_aligned_reads.bam -sorted > gene_coverage.txt

coverageBed -a chr21_22_annotation_sorted.bed -b sorted_aligned_reads_only.bam > gene_coverage.txt

# Task5_coverage_noramlisation.py gets the number of reads mapped per exon per gene
# and computes the estimate number of reads per gene.
# The final output file (reads_per_gene_normalisations.txt) will have RPM and RPKM 
# normalisations in the last two columns

chmod +x Task5_coverage_normalisation.py
./Task5_coverage_normalisation.py gene_coverage.txt

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
