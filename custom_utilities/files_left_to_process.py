all_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\recovery\\all.txt"
done_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\recovery\\done.txt"
left_file = "C:\\Users\\jorsm\\OneDrive\\masters_thesis\\VCFs\\vcfs\\recovery\\left.txt"

left=[]

def file_list_parser(file):
    items=[]
    with open(file, 'r') as list:
        for line in list:
            line = line.strip('\n')
            items.append(line)
    return items

all_items = file_list_parser(all_file)
done_items = file_list_parser(done_file)
print(f"all:",all_items)
print(f"done:",done_items)

for item in all_items:
    if item not in done_items:
        left.append(item)

print(f"left:",left)

with open(left_file, 'w') as left_file_writer:
    for item in left:
        #item = item.strip("vcf_")
        item = item + '\n'
        left_file_writer.write(item)
