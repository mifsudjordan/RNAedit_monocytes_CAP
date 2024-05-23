#!/usr/bin/env python3

input_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\paired_variants.txt"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\paired_samples_variants.csv"

ids = []
individuals = []
ind_acute_only = []
acute_var_id = []
rec_var_id = []
acute_only_var = []
rec_only_var = []
percentage = []
perc_acute_only = []

with open(input_file, "r") as variants:
    # Filling up the list of ids
    for line in variants:
        field = line.split("\t")
        if field[7].strip('\n') not in ids:
            ids.append(field[7].strip('\n'))
    variants.seek(0)
    print(f"There are ", len(ids)," individuals to go through")

    # Filling lists of acute and recovery variants and their
    # percentage occurance in whole group and ids
    for id in ids:
        for line in variants:
            field = line.split("\t")
            if field[7].strip('\n') == id:
                if field[5].strip('\t') == "Hospital admission":
                    acute_var_id.append(field[0].strip('\t'))
                    percentage.append(field[4].strip('\t'))
                    individuals.append(field[7].strip('\n'))
                else:
                    rec_var_id.append(field[0].strip('\t'))

        variants.seek(0)

        # Finding variants unique to acute phase in the same individual
        counter = 0
        for variant in acute_var_id:
            if variant not in rec_var_id:
                acute_only_var.append(variant)
                perc_acute_only.append(percentage[counter])
                ind_acute_only.append(individuals[counter])
            counter = counter + 1

        # Finding variants unique to recovery phase in the same individual
        for variant in rec_var_id:
            if variant not in acute_var_id:
                rec_only_var.append(variant)

        # Initialising lists
        acute_var_id = []
        rec_var_id = []
        percentage = []

# Initialising lists for the next code-block
unique_acute_only_var = []
counts = []
counts_2 = []


# Count of occurances of each acute-only-variant in the whole group
for variant in acute_only_var:
    counts.append(acute_only_var.count(variant))

# Remove duplicates from list of variants
counter = 0
for variant in acute_only_var:
    if variant not in unique_acute_only_var:
        unique_acute_only_var.append(variant)
        counts_2.append(counts[counter])
    counter = counter + 1

print(f"Total acute-only variants:", len(unique_acute_only_var))

# Combine the two lists into a list of tuples, and sort them based on counts_2 in descending order
sorted_data = sorted(zip(unique_acute_only_var, counts_2), key=lambda x: x[1], reverse=True)

# Unpack the sorted data into separate lists
sorted_unique_acute_only_var, sorted_counts_2 = zip(*sorted_data)

# Sorted_unique_acute_only_var contains the sorted strings,
# and sorted_counts_2 contains the sorted integers.

# Writing output file
with open(input_file, "r") as records_file:
    list_of_ids = ''
    with open(output_file, "w") as output:
        header = "unique_var_id\tcount_in_acute_only\tcount_in_rec_only\tfound_in_id\n"
        output.write(header)

        counter = 0
        total_counts = 0
        for record in sorted_unique_acute_only_var:
            for line in records_file:
                field = line.split('\t')
                if field[0] == record:
                    list_of_ids = list_of_ids + field[7].strip('\n') + ','
                    total_counts = total_counts + 1
            records_file.seek(0)
            list_of_ids = list_of_ids[:-1] # Removes last comma from string
            counts_in_acute = str(sorted_counts_2[counter])
            counts_in_rec = str(total_counts - sorted_counts_2[counter])
            line_to_write = record + "\t" + counts_in_acute + "\t" + counts_in_rec + '\t' + list_of_ids + "\n"
            output.write(line_to_write)
            counter = counter + 1
            list_of_ids = ''
            total_counts = 0
#print(sorted_unique_acute_only_var)
#print(sorted_counts_2)
