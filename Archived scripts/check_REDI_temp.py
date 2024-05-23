#!/usr/bin/env python3

import sys
vcf_file = sys.argv[1]

red_database_directory ="/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/indexedAtoI/"
red_database_file = "/mnt/c/Users/jorsm/OneDrive/masters_thesis/VCFs/REDIDatabase_hg38.txt"
#vcf_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\hosp_admis\\vcf_SRR12926209\\SRR12926209_AtoI.vcf"


chromosome = []
positions = []

chromosome_vcf = []
positions_vcf = []
records_vcf = []
records_database = []

prefix = 'chr'
matches = []
n = 0

with open(vcf_file, "r") as vcf:
    for line in vcf:
        if not line.startswith("#"):
            field = line.split('\t')
            if len(field[0]) < 3:
                chr_name = prefix + field[0]
                if chr_name not in chromosome_vcf:
                    chromosome_vcf.append(chr_name)


match_counter = 0
total_matches= 0
for vcf_chromosome in chromosome_vcf:
    check_chr = vcf_chromosome.strip(prefix)
    print(f"Matches in Chr ",check_chr)
    coordinates = []
    coordinates_vcf = []
    if check_chr == "MT":
        vcf_chromosome = "chrM"
    path_to_red_database_index = red_database_directory + "AtoIedits_in_" + vcf_chromosome + ".txt"
    with open(path_to_red_database_index, "r") as index:
        for line in index:
            coordinates.append(line.strip("\n"))
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if not line.startswith("#"):
                field = line.split('\t')
                if field[0] == check_chr:
                    coordinates_vcf.append(field[1])
    for record in coordinates_vcf:
        if record in coordinates:
            match_counter = match_counter + 1
            #print(record)
    print(match_counter)
    total_matches = total_matches + match_counter
    match_counter = 0
print("Total matches to REDI database")
print(total_matches)