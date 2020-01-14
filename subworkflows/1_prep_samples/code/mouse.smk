
checkpoint download_mouse:
    input:
        list="data/mouse/SRR_Acc_List.txt",
        sh="code/download.sh"
    output:
        dir=dir("data/mouse/raw")
    params:
        tar="data/mouse/raw/StabilityNoMetaG.tar",
        url="http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar"
    benchmark:
        "benchmarks/mouse/download.txt"
    shell:
        """
        wget -N -P {output.dir} {params.url}
        tar -xvf {params.tar} -C {output.dir}
        rm {params.tar}
        """
        
def get_mouse_fastq(wildcards):
    dir = checkpoints.download_mouse.get().output.dir
    filenames, = glob_wildcards(os.path.join(dir, "{filename}"))
    return expand("{dir}/{filename}", dir=dir, filename=filenames)

rule names_file_mouse:
    input:
        R="code/mouse.R",
        files=get_mouse_fastq
    output:
        file="data/mouse/mouse.files"
    benchmark:
        "benchmarks/mouse/names_file.txt"
    script:
        "{input.R}"