#!/usr/bin/env python3

sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers.txt"
occursinlist = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\occuts_in_interesting_variants.txt"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\interesting_variants_samples.txt"
groups ='ar'
with open (occursinlist, "r") as occursinfile:
    for line in occursinfile:
        occurs_in = line.strip('\n')
        srr_numbers = []
        with open(sample_list, "r") as samples:
            for line in samples:
                sample = line.strip('\n')
                srr_numbers.append(sample)
        matching_samples = []
        counter = 0
        for digit in occurs_in:
            if digit in groups:
                matching_samples.append(srr_numbers[counter])
            counter = counter + 1
        with open(output_file, 'a') as output:
            line_to_write = str(matching_samples) + '\n'
            output.write(line_to_write)



