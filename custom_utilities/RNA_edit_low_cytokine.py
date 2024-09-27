from scipy.stats import mannwhitneyu
import numpy as np
#from statsmodels.stats.multitest import multipletests

cytokine_data = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\Cytokine_release_data.txt"
sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers.txt"
acute_sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers_hosp.txt"
variant_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\ctou_ann_stats_strict_p.csv"
output_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\LPS_IL_1b_ctou_acute_strict_pvalues.txt"

occurs_in = ''
groups ='a'
p_values = []
mean_values = []
median_values = []
mean_il6_values = []
median_il6_values = []

# Getting list of SRR numbers
srr_numbers = []
with open(sample_list, "r") as samples:
    for line in samples:
        sample = line.strip('\n')
        srr_numbers.append(sample)

# Getting list of acute SRR numbers
srr_acute = []
with open(acute_sample_list, "r") as acute_samples:
    for line in acute_samples:
        sample = line.strip('\n')
        srr_acute.append(sample)

with open(variant_list, "r") as variant_file:
        for line in variant_file:
            field = line.split('\t')
            if field[0] != "Chr":
                occurs_in = field[10]

                # Converting occurs_in into a list of acute SRR numbers with that variant
                matching_samples = []
                counter = 0
                for digit in occurs_in:
                    if digit in groups:
                        matching_samples.append(srr_numbers[counter])
                    counter = counter + 1

                # List of acute SRR numbers without that variant
                non_matching_samples = []
                counter = 0
                for srr in srr_acute:
                    if srr not in matching_samples:
                        non_matching_samples.append(srr)

                match_il10_values = []
                non_match_il10_values = []
                with open(cytokine_data, "r") as cytokines:
                    for line in cytokines:
                        if line.startswith("SRR"):
                            field = line.split('\t')
                            if field[0] in matching_samples:
                                match_il10_values.append(field[6])
                            if field[0] in non_matching_samples:
                                non_match_il10_values.append(field[6])

                # cleaning lists from empty values and converting to a list of float numbers
                temp_list = []
                temp_list = [float(x) for x in match_il10_values if x != '']
                match_il10_values = []
                match_il10_values = temp_list

                temp_list = []
                temp_list = [float(x) for x in non_match_il10_values if x != '']
                non_match_il10_values = []
                non_match_il10_values = temp_list


                # Assuming match_il10_values and non_match_il10_values are your two lists
                if len(match_il10_values) > 1:
                    match_mean = np.mean(match_il10_values)
                    match_median = np.median(match_il10_values)
                    mean_values.append(match_mean)
                    median_values.append(match_median)

                    if len(non_matching_samples) > 1:
                        stat, p_value = mannwhitneyu(match_il10_values, non_match_il10_values, alternative='two-sided')
                        p_values.append(p_value)
                    else:
                        p_values.append("NA")

                else:
                    p_values.append("NA")
                    mean_values.append("NA")
                    median_values.append("NA")

                #print(f"Statistic: {stat}")
                #print(f"P-value: {p_value}")

#p_values = np.array(p_values)

# Perform Benjamini-Hochberg adjustment
#_, bh_adjusted_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')


#print("BH adjusted p-values:", bh_adjusted_p_values)

counter = 0
with open(output_file, "w") as output:
    for value in p_values:

        linetowrite = ''
        item = str(value)
        mean = str(mean_values[counter])
        median = str(median_values[counter])
        linetowrite = item + '\t' + mean + '\t' + median + '\n'
        output.write(linetowrite)
        counter +=1