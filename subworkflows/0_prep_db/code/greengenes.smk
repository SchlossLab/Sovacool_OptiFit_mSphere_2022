subworkflow silva:
    snakefile:
        "code/silva.smk"

mothur = "mothur '#set.dir(input=data/gg/, output=data/gg/); set.logfile(name={log}); "

rule download_gg:
    output:
        fasta="data/gg/gg_13_8_99.fasta",
        tax="data/gg/gg_13_8_99.tax"
    shell:
        """
        wget -N -P data/gg/ http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
        wget -N -P data/gg/ http://www.mothur.org/w/images/1/19/Gg_13_8_99.refalign.tgz
        tar -xzvf data/gg/Gg_13_8_99.taxonomy.tgz -C data/gg/
        tar -xzvf data/gg/Gg_13_8_99.refalign.tgz -C data/gg/
        """

rule get_gg_bact:
    input:
        fasta=rules.download_gg.output.fasta,
        tax=rules.download_gg.output.tax
    output:
        fasta="data/gg/gg.bacteria.fasta",
        tax="data/gg/gg.bacteria.tax",
    shell:
        """
        {mothur}
        get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria)'
        mv data/gg/gg_13_8_99.pick.fasta data/gg/gg.bacteria.fasta
        mv data/gg/gg_13_8_99.pick.tax data/gg/gg.bacteria.tax
        """

rule align_gg:
    input:
        silva_ref=silva("data/silva/silva.seed_v132.align"),
        fasta=rules.get_gg_bact.output.fasta
    output:
        align="data/gg/gg.bacteria.align"
    shell:
        """
        {mothur}
        align.seqs(candidate={input.fasta}, template={input.silva_ref})
        '
        """

rule get_gg_full:
    input:
        align=rules.align_gg.output.align
    output:
        screen="data/gg/gg.bacteria.good.align",
        filter="data/gg/gg.bacteria.good.filter.fasta",
        names="data/gg/gg.bacteria.good.filter.names",
        unique="data/gg/gg.bacteria.good.filter.unique.fasta",
        preclust_fasta="data/gg/gg.bacteria.good.filter.unique.precluster.fasta",
        preclust_names="data/gg/gg.bacteria.good.filter.unique.precluster.names",
    shell:
        """
        {mothur}
        screen.seqs(fasta={input.align}, start=2000, end=414788)
        filter.seqs(fasta={output.screen}, trump=., vertical=T)
        unique.seqs(fasta={output.filter})
        pre.cluster(fasta={output.unique}, name={output.names}, diffs=2)
        '
        """
