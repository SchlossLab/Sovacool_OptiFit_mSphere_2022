rule download_marine:
    input:
        list="data/marine/SRR_Acc_List.txt",
        sh="code/download.sh"
    output:
        fastq=fastq_filenames["marine"]
    benchmark:
        "benchmarks/marine/download.txt"
    shell:
        "{input.sh} {input.list} data/marine/raw/"

rule names_file_marine:
    input:
        R="code/marine.R",
        files=rules.download_marine.output.fastq
    output:
        file="data/marine/marine.files"
    benchmark:
        "benchmarks/marine/names_file.txt"
    script:
        "marine.R"