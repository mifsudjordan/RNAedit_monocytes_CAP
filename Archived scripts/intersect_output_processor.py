#!/usr/bin/env python3

sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers.txt"
intersect_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\all_ctou_vcf_compressed\\intersect_ctou_all.vcf"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\all_ctou_vcf_compressed\\modified_intersect_ctou_all.csv"

srr_numbers = []
with open(sample_list, "r") as samples:
    for line in samples:
        sample = line.strip('\n')
        srr_numbers.append(sample)

srr_list = ''
for srr in srr_numbers:
    srr_list = srr_list + srr + '\t'

srr_list = srr_list[:-1]

with open(output_file, "w") as output:
    header = "unique_var_id" + "\t" + srr_list + "\n"
    output.write(header)

    with open(intersect_file, "r") as intersects:
        for line in intersects:
            field = line.split('\t')
            chr = field[0]
            pos = field[1]
            ref = field[2]
            alt = field[3]
            found_in = field[4].strip('\n')
            modified_found_in = ''
            for digit in found_in:
                if digit == "1":
                    modified_found_in = modified_found_in + "CtoU" + "\t"
                if digit == "0":
                    modified_found_in = modified_found_in + "none" + "\t"

            line_to_write = chr + ":" + pos + ref + alt + "\t" + modified_found_in + "\n"

            output.write(line_to_write)




