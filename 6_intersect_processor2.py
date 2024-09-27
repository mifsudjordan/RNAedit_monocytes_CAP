#!/usr/bin/env python3

time_points = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\timepoints.txt"
intersect_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\annovar_output\\intersect_atoi_strict_ann.hg38_multianno.txt"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\annovar_output\\atoi_ann_stats_strict.csv"


# Getting list of time_points as a list of symbols: a, r, and r.
group = []
with open(time_points, "r") as time_point_list:
    for line in time_point_list:
        timepoint = line.strip('\n')
        group.append(timepoint[0])

total_a = group.count("a")
total_c = group.count("c")
total_r = group.count("r")
print(total_a, total_r, total_c)

# The following code takes every annotated variant / RNA edit and
# counts the occurrences of each variant per time-point (recovery, acute or controls)
# calculates the percentage occurrence in each group.

counter = 0
additional_columns = "perc_total\ta_count\tc_count\tr_count\ta_noncount\tc_noncount\tr_noncount\tperc_acute\tperc_con\tperc_rec\n"
with open(output_file, "w") as output:
    with open(intersect_file, "r") as intersects:
        for line in intersects:
            if line.startswith("C"):
                header = line.strip('\n') + '\t' + additional_columns
                output.write(header)
            else:
                field = line.split('\t')
                found_in = field[10].strip('\n')
                matches = found_in.count("1")
                total_samples = len(found_in)
                perc_total = round(matches / total_samples * 100, 2)

                #Coverting 1s to symbols of group
                occurs_in = '' #empty string for list of occurances: "a", "r" and "c" and "0".
                for sample in found_in:
                    if sample == ("1"):
                        occurs_in = occurs_in + str(group[counter])
                    else:
                        occurs_in = occurs_in + "0"
                    counter = counter + 1
                counter = 0 #restarting counter
                occurs_in = occurs_in + '\t'

                #Updating record with occurs_in and perc_total
                record = field[:-1]
                line_to_write = '' #emptying line to write
                for item in record:
                    line_to_write = line_to_write + str(item) + "\t"
                line_to_write = line_to_write + occurs_in + str(perc_total) + '\t'

                #Calculating stats for acute, controls and recovery
                a_matches = occurs_in.count("a")
                c_matches = occurs_in.count("c")
                r_matches = occurs_in.count("r")
                a_non_matches = str(total_a - a_matches) + '\t'
                c_non_matches = str(total_c - c_matches) + '\t'
                r_non_matches = str(total_r - r_matches) + '\t'

                matches_string = str(a_matches) + '\t' + str(c_matches) + '\t' + str(r_matches) + '\t'
                non_matches_string = a_non_matches + c_non_matches + r_non_matches

                line_to_write = line_to_write + matches_string + non_matches_string

                perc_a = str(round(a_matches / total_a * 100, 2)) + '\t'
                perc_c = str(round(c_matches / total_c * 100, 2)) + '\t'
                perc_r = str(round(r_matches / total_r * 100, 2))

                #Completing line to write
                line_to_write = line_to_write + perc_a + perc_c + perc_r + '\n'
                output.write(line_to_write)



