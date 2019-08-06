rule download_rdp:
    output:
        tar=temp("data/rdp/Trainset16_022016.rdp.tgz"),
        fasta="data/rdp/rdp.fasta",
        tax="data/rdp/rdp.tax"
    params:
        fasta="data/rdp/trainset16_022016.rdp/trainset16_022016.rdp.fasta",
        tax="data/rdp/trainset16_022016.rdp/trainset16_022016.rdp.tax"
    shell:
        """
        wget -N -P data/rdp/ https://www.mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz
        tar xvzf {output.tar} -C data/rdp/
        mv {params.fasta} {output.fasta}
        mv {params.tax} {output.tax}
        rm -rf data/rdp/trainset16_022016.rdp/
        """
