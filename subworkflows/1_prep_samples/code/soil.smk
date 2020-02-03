rule names_file_soil:
    input:
        R="code/soil.R",
        files=expand("data/soil/raw/{SRA}_{i}.fastq.gz", SRA=sra_list["soil"], i=(1,2))
    output:
        file="data/soil/soil.files"
    benchmark:
        "benchmarks/soil/names_file.txt"
    shell:
        "Rscript {input.R}"