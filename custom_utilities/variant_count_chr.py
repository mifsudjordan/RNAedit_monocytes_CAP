#!/usr/bin/env python3

import os

vcf_file_path = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\all_atoi_strict\\"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\all_atoi_strict\\AtoI_strict_count_chr.csv"

# creating the header with all the chromosomes
chr = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
       "11", "12", "13", "14", "15", "16", "17", "18",
       "19", "20", "21", "22", "X", "Y", "MT"]
header_prefix = "AtoI_Chr"
header = ''
for item in chr:
    header = header + header_prefix + item + "\t"
header = header[:-1]
header = header + "\n"
header = "Run" + "\t" + "total_atoi" + "\t" + header
print(header)

with open(output_file, "w") as output:
    output.write(header)

file_suffix = "_AtoI_strict.vcf"

for file in os.listdir(vcf_file_path):
    if file.endswith('.vcf'):
        chromosomes = []
        counts = []
        variant_counter = 0
        filename_fields = file.split("_")
        srr = filename_fields[0]
        full_file_path = vcf_file_path + srr + file_suffix

        #Parsing of the vcf file to count variants per chromosome
        with open(full_file_path, "r") as file:
            for line in file:
                if not line.startswith("#"):
                    field = line.split("\t")
                    if field[0] not in chromosomes:
                        if chromosomes: # checks if chromosomes list is still empty
                            counts.append(variant_counter) # this is executed only after parsing through first chromosome
                        variant_counter = 1
                        chromosomes.append(field[0])
                    else:
                        variant_counter = variant_counter + 1
            counts.append(variant_counter)

        # re-ordering the chromosome in ascending order and adding 0 whenever a chromosome has no variant
        ordered_counts = []
        for chromosome in chr:
            if chromosome in chromosomes:
                ordered_counts.append(counts[chromosomes.index(chromosome)])
            else:
                ordered_counts.append(0)

        total_ctou = str(sum(counts))
        # preparing line to write
        line_to_write = ""
        for i in ordered_counts:
            line_to_write = line_to_write + str(i) + "\t"
        line_to_write = line_to_write[:-1]
        line_to_write = line_to_write + "\n"
        line_to_write = srr + "\t" + total_ctou + "\t" + line_to_write

        with open(output_file, "a") as output:
            output.write(line_to_write)
