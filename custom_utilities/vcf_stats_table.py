#!/usr/bin/env python3

stats_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\controls\\count_CtoUchr_only2.txt"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\controls\\stats_table_controls3.txt"

srr = []
headers = []
quantity = []
counter = 0
record = 0
line_counter = 0
line_counter_actual = 0
header_written = False
with open(stats_file, "r") as file:
    for line in file:
        line_counter = line_counter + 1

with open(stats_file, "r") as file:
    for line in file:
        line_counter_actual = line_counter_actual + 1
        line = line.strip('\n')
        if counter == 2:
            line_to_write = srr[record] + '\t'
            quantities = len(quantity)
            for n in range(quantities):
                line_to_write = line_to_write + quantity[n] + '\t'
            line_to_write = line_to_write + '\n'
            if not header_written:
                headers_to_write = 'SRR number\t'
                headers_count = len(headers)
                for n in range(headers_count):
                    headers_to_write = headers_to_write + headers[n] + '\t'
                headers_to_write = headers_to_write + '\n'
                with open(output_file, "a") as table_file:
                    table_file.write(headers_to_write)
                header_written = True
            with open(output_file, "a") as table_file:
                table_file.write(line_to_write)
                counter = 1
                record = record + 1
                quantity = []
                line_to_write = ''

        if line.startswith("SRR"):
            srr.append(line)
            counter = counter + 1
        if not line.startswith("SRR"):
            if not line.isnumeric():
                headers.append(line)
            else:
                quantity.append(line)

        # ---------------to write the last record
        if line_counter_actual == line_counter:
            line_to_write = srr[record] + '\t'
            quantities = len(quantity)
            for n in range(quantities):
                line_to_write = line_to_write + quantity[n] + '\t'
            line_to_write = line_to_write + '\n'
            if not header_written:
                headers_to_write = 'SRR number\t'
                headers_count = len(headers)
                for n in range(headers_count):
                    headers_to_write = headers_to_write + headers[n] + '\t'
                headers_to_write = headers_to_write + '\n'
                with open(output_file, "a") as table_file:
                    table_file.write(headers_to_write)
                header_written = True
            with open(output_file, "a") as table_file:
                table_file.write(line_to_write)
    no_of_srr = len(srr)
print(srr)
print(headers)
print(quantity)






