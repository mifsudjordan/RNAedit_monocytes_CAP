#!/usr/bin/env python3

import sys

red_database_directory ="C:/Users/jorsm/OneDrive/masters_thesis/VCFs/indexedAtoI/"
redi_database_file = "C:/Users/jorsm/OneDrive/masters_thesis/VCFs/REDIDatabase_hg38.txt"
intersect_file = "C:/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/all_atoi_strict_compressed/intersect_atoi_strict.vcf"

chr = ''
unique_var_id = ''
my_edits = []
common_indexed_filename = 'AtoIedits_in_chr'
unique_edit_id = ''
matches = 0
checked = 0

chromosomes = []

with open(intersect_file, "r") as int_file:
    for line in int_file:
        field = line.split('\t')
        chr = field[0]

        # changes naming convention of the mitochondrial chromosome
        if chr == "MT":
            chr = "M"

        if chr not in chromosomes:
            chromosomes.append(chr)

        unique_var_id = chr + "." + field[1]
        my_edits.append(unique_var_id)
print(my_edits)
for i in chromosomes:
    indexed_file = red_database_directory + common_indexed_filename + i + '.txt'
    print(indexed_file)
    print(i)
    with open(indexed_file, "r") as indexed_db:
        redi_edits = []
        for line in indexed_db:
            redi_edits.append(line.strip("\n"))
        for edit in my_edits:
            potential_edit = edit.split(".")
            if potential_edit[0] == i:
                if potential_edit[1] in redi_edits:
                    matches +=1
            sys.stdout.write(f"\rMatches: {matches}")
            sys.stdout.flush()
print(matches)