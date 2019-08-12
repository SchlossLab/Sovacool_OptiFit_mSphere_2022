rule download_rdp:
    output:
        tar=temp("data/rdp/Trainset16_022016.pds.tgz"),
        fasta="data/rdp/rdp.fasta",
        tax="data/rdp/rdp.tax"
    params:
        fasta="data/rdp/trainset16_022016.pds/trainset16_022016.pds.fasta",
        tax="data/rdp/trainset16_022016.pds/trainset16_022016.pds.tax"
    shell:
        """
        wget -N -P data/rdp/ https://mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz
        tar xvzf {output.tar} -C data/rdp/
        mv {params.fasta} {output.fasta}
        mv {params.tax} {output.tax}
        rm -rf data/rdp/trainset16_022016.pds/
        """
