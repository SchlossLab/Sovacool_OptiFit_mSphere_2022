
with open("data/mouse/SRR_Acc_List.txt", 'r') as file:
    mouse_filenames = [f"data/mouse/raw/{line.strip()}" for line in file]

rule download_mouse:
    input:
        #list="data/mouse/SRR_Acc_List.txt",
        sh="code/download.sh"
    output:
        files=mouse_filenames
    params:
        tar="data/mouse/StabilityNoMetaG.tar",
        url="http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar"
    benchmark:
        "benchmarks/mouse/download.txt"
    shell:
        """
        wget -N -P data/mouse/ {params.url}
        tar -xvf {params.tar} -C data/mouse/raw/
        rm {params.tar}
        """

rule names_file_mouse:
    input:
        files=rules.download_mouse.output.files
    output:
        file="data/mouse/mouse.files"
    benchmark:
        "benchmarks/mouse/names_file.txt"
    run:
        input_files = snakemake.input['files']
        output_file = snakemake.output['file']

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