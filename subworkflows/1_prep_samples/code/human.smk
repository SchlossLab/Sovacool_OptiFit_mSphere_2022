rule download_human:
    input:
        list="data/human/SRR_Acc_List.txt",
        sh="code/download.sh"
    output:
        fastq=fastq_filenames["human"]
    benchmark:
        "benchmarks/human/download.txt"
    shell:
        "{input.sh} {input.list} data/human/raw/"

rule names_file_human:
    input:
        R="code/human.R",
        files=rules.download_human.output.fastq
    output:
        file="data/human/human.files"
    benchmark:
        "benchmarks/human/names_file.txt"
    script:
        "human.R"