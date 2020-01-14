rule download_soil:
    input:
        list="data/soil/SRR_Acc_List.txt",
        sh="code/download.sh"
    output:
        fastq=fastq_filenames["soil"]
    benchmark:
        "benchmarks/soil/download.txt"
    shell:
        "{input.sh} {input.list} data/soil/raw/"

rule names_file_soil:
    input:
        R="code/soil.R",
        files=rules.download_soil.output.fastq
    output:
        file="data/soil/soil.files"
    benchmark:
        "benchmarks/soil/names_file.txt"
    script:
        "{input.R}"