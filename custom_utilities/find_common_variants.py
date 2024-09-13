#!/usr/bin/env python3

sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers_controls.txt"
intersect_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\controls\\intersect\\intersect_controls.vcf"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\controls\\intersect\\common_variants_controls.txt"

srr_numbers = []
with open(sample_list, "r") as samples:
    for line in samples:
        sample = line.strip('\n')
        srr_numbers.append(sample)

# This code takes the intersect output and finds the percentage occurrence
# of each variant.

with open(output_file, "w") as output:
    header1 = "#Processed file: " + intersect_file + '\n'
    header2 = "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "matches" + "\t" + "percentage" + "\t" + "Found in" + "\n"
    output.write(header1)
    output.write(header2)

    with open(intersect_file, "r") as intersects:
        for line in intersects:
            matching_samples = []
            field = line.split('\t')
            chr = field[0]
            pos = field[1]
            ref = field[2]
            alt = field[3]
            found_in = field[4].strip('\n')
            matches = found_in.count("1")
            total_samples = len(found_in)
            percentage = round(matches / total_samples * 100, 2)

            counter = 0
            for digit in found_in:
                if digit == "1":
                    matching_samples.append(srr_numbers[counter])
                counter = counter + 1

            matching_string = ''
            for srr in matching_samples:
                matching_string = matching_string + srr + ","

            line_to_write = chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + str(matches) + "\t" + str(percentage) + "\t" + matching_string + "\n"

            output.write(line_to_write)




