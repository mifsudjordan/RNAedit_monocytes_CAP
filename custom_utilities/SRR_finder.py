cytokine_data = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\Cytokine_release_data.txt"
sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers.txt"

occurs_in = "a000000000000000a00000000000000a0a00000000000000000000000000a000000000r000000000a00000000000000000000000000000000000000000000000000000000a0000000a000a00a000000a00000r"
groups ='ar'

# Getting list of SRR numbers
srr_numbers = []
with open(sample_list, "r") as samples:
    for line in samples:
        sample = line.strip('\n')
        srr_numbers.append(sample)

# Converting occurs_in into a list of acute SRR numbers with that variant
matching_samples = []
counter = 0
for digit in occurs_in:
    if digit in groups:
        matching_samples.append(srr_numbers[counter])
    counter = counter + 1
print(matching_samples)




