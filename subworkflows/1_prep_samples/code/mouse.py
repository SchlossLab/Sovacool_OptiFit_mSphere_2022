import collections
import os

if 'snakemake' in locals() or 'snakemake' in globals():
    print("Using filepaths from snakemake")
    input_files = snakemake.input['files']
    output_file = snakemake.output['file']
else:
    print("Using default filepaths")
    input_files = os.listdir("data/mouse/raw")
    output_file = 'data/mouse/mouse.files'

print("Mapping names to files")
names = collections.defaultdict(list)
for file in sorted(input_files):
    name = file.split('/')[-1].split('_')[0]
    if "mock" not in name.lower():
        names[name].append(file)

print(f"Writing {output_file}")
with open(output_file, 'w') as file:
    for name, files in names.items():
        file_list = '\t'.join(files)
        file.write(f"{name}\t{file_list}\n")