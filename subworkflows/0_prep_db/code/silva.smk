rule download_silva:
    output:
        tar="data/silva/Silva.seed_v132.tgz",
        fasta="data/silva/silva.fasta",
        tax="data/silva/silva.tax"
    params:
        fasta="data/silva/silva.seed_v132.align",
        tax="data/silva/silva.seed_v132.tax"
    shell:
        """
        wget -N -P data/silva/ http://www.mothur.org/w/images/7/71/Silva.seed_v132.tgz
        tar xzvf {output.tar} -C data/silva/
        mv {params.fasta} {output.fasta}
        mv {params.tax} {output.tax}
        """

rule rename_silva_bact:
    input:
        "data/silva/silva.bacteria.fasta"
    output:
        "data/silva/silva.bacteria.align"
    shell:
        "cp {input} {output}"
