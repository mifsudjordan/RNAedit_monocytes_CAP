#!/usr/bin/env python3

stats_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\recovery\\variant_stats_recovery_strict.txt"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\recovery\\variant_stats_recovery_strict.csv"

headers_to_write = "Run\tuncommon_snps_qc\tAtoIandCtoU\tAtoI\tCtoU\tAtoI_Ch1matches\tAtoI_Ch10matches\t \\" \
                   "AtoI_Ch11matches\tAtoI_Ch12matches\tAtoI_Ch13matches\t \\" \
                   "AtoI_Ch14matches\tAtoI_Ch15matches\tAtoI_Ch16matches\tAtoI_Ch17matches\tAtoI_Ch18matches\t \\" \
                   "AtoI_Ch19matches\tAtoI_Ch2matches\tAtoI_Ch20matches\tAtoI_Ch21matches\tAtoI_Ch22matches\t \\" \
                   "AtoI_Ch3matches\tAtoI_Ch4matches\tAtoI_Ch5matches\tAtoI_Ch6matches\tAtoI_Ch7matches\t \\" \
                   "AtoI_Ch8matches\tAtoI_Ch9matches\tAtoI_ChMTmatches\tAtoI_ChXmatches\tAtoI_ChYmatches\t \\" \
                   "totalmatches\n"

#both_symbol = "||"
#ctou_symbol = "No. of C to U:"
quantity = []
quantity_counter = 0
line_to_write = ''
counter = 0

# This code takes the output of "script_variant_filtering_counter_v2".sh
# containing all the SRR numbers of a time-point and
# all it's RNA edit counts per chromosome
# and converts it into a tab-delimited csv file.

with open(stats_file, "r") as file:
    with open(output_file, "a") as table_file:
        table_file.write(headers_to_write)
        for line in file:
            line = line.strip("\n")
            if line.startswith("SRR"):
                run = line

            if line.startswith("Total"): #conditional in case there is no record of Y chromosome
                if quantity_counter == 28:
                    quantity_counter = quantity_counter + 1
                    quantity.append("0")

            if line.isnumeric():
                quantity.append(line)
                quantity_counter = quantity_counter + 1

            if quantity_counter == 30:
                for item in quantity:
                    line_to_write = line_to_write + item + '\t'
                line_to_write = line_to_write[:-1]
                line_to_write = line_to_write + '\n'
                line_to_write = run + '\t' + line_to_write
                quantity_counter = 0
                table_file.write(line_to_write)
                line_to_write = ''
                quantity = []
                counter = counter + 1
print(f"No. of samples: ", counter)


