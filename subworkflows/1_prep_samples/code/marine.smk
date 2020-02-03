rule names_file_marine:
    input:
        R="code/marine.R",
        files=expand("data/marine/raw/{SRA}_{i}.fastq.gz", SRA=sra_list["marine"], i=(1,2))
    output:
        file="data/marine/marine.files"
    benchmark:
        "benchmarks/marine/names_file.txt"
    script:
        "code/marine.R"