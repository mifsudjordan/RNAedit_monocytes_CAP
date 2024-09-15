#!/usr/bin/env python3

#import sys
#intersect_file = sys.argv[1]

intersect_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\all_atoi_strict_compressed\\intersect_atoi_strict.vcf"

#output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\all_ctou_strict_compressed\\intersect_annovar_ctou_strict.vcf"
output_file = "/Users/jorsm/OneDrive/masters_thesis/VCFs/vcfs/all_atoi_strict_compressed/intersect_annovar_atoi_strict.vcf"

with open(output_file, "w") as output:
    with open(intersect_file, "r") as intersect:
        for line in intersect:
            field = line.split('\t')
            chr = field[0] + '\t'
            pos = field[1] + '\t'
            pos2 = field[1] + '\t'
            ref = field[2] + '\t'
            alt = field[3] + '\t'
            found_in = field[4]

            line_to_write = chr + pos + pos2+ ref + alt + found_in

            output.write(line_to_write)