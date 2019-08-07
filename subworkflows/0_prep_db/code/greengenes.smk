rule download_gg:
    output:
        fasta="data/gg/gg.fasta",
        tax="data/gg/gg.tax"
    params:
        fasta="data/gg/gg_13_8_99.fasta",
        tax="data/gg/gg_13_8_99.gg.tax"
    shell:
        """
        wget -N -P data/gg/ http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
        wget -N -P data/gg/ http://www.mothur.org/w/images/1/19/Gg_13_8_99.refalign.tgz
        tar -xzvf data/gg/Gg_13_8_99.taxonomy.tgz -C data/gg/
        tar -xzvf data/gg/Gg_13_8_99.refalign.tgz -C data/gg/
        mv {params.fasta} {output.fasta}
        mv {params.tax} {output.tax}
        """
