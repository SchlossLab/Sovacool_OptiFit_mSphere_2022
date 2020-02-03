rule names_file_human:
    input:
        R="code/human.R",
        files=expand("data/human/raw/{SRA}_{i}.fastq.gz", SRA=sra_list["human"], i=(1,2))
    output:
        file="data/human/human.files"
    benchmark:
        "benchmarks/human/names_file.txt"
    script:
        "code/human.R"