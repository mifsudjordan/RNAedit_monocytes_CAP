cytokine_data = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\masters-analysis\\Cytokine_release_data.txt"
sample_list = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\srr_numbers.txt"

occurs_in = "a0000a000ra00r00a0000r000000a0000a00r00000ara0r0000a000r0r000r0000r000r000000r00a0000a0a00000000000c000c000000000000000000000000c0c000000a0000ra0a00000r0000r00a0a0a0r"
groups ='arc'

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
print(len(matching_samples))




