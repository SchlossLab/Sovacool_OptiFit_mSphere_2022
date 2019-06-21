subworkflow silva:
    snakefile:
        "code/silva.smk"

mothur = "mothur '#set.dir(input=data/rdp/, output=data/rdp/); set.logfile(name={log}) "

rule download_rdp:
    output:
        tar=temp("data/rdp/Trainset16_022016.rdp.tgz"),
        fasta="data/rdp/rdp.fasta",
        tax="data/rdp/rdp.tax"
    shell:
        """
        wget -N -P data/rdp/ https://www.mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz
        tar xvzf {output.tar} -C data/rdp/
        mv data/rdp/trainset16_022016.rdp/trainset16_022016.rdp.fasta {output.fasta}
        mv data/rdp/trainset16_022016.rdp/trainset16_022016.rdp.tax {output.tax}
        rm -rf data/rdp/trainset16_022016.rdp/
        """

rule get_rdp_bact:
    input:
        fasta=rules.download_rdp.output.fasta,
        tax=rules.download_rdp.output.fasta
    output:
        fasta="data/rdp/rdp.bacteria.fasta",
        tax="data/rdp/rdp.bacteria.tax"
    log:
        "logfiles/rdp/get_rdp_bact.log"
    shell:
        """
        {mothur}
        get.lineage(fasta={input.fasta}, tax={input.tax}, taxon=Bacteria)'
        mv data/rdp/rdp.pick.fasta {output.fasta}
        mv data/rdp/rdp.pick.tax {output.tax}
        """

rule get_full_length_rdp:
    input:
        fasta=rules.get_rdp_bact.output.fasta,
        tax=rules.get_rdp_bact.output.tax,
        silva_ref=silva("data/silva/silva.seed_v132.align")
    output:
        align="data/rdp/rdp.bacteria.align",
        report="data/rdp/rdp.bacteria.report",
        screen="data/rdp/rdp.bacteria.good.align",
        filter="data/rdp/rdp.bacteria.good.filter.fasta",
        names="data/rdp/rdp.bacteria.good.filter.names",
        unique="data/rdp/rdp.bacteria.good.filter.unique.fasta",
        preclust_fasta="data/rdp/rdp.bacteria.good.filter.unique.precluster.fasta",
        preclust_names="data/rdp/rdp.bacteria.good.filter.unique.precluster.names"
    shell:
        """
        {mothur}
        align.seqs(candidate={input.fasta}, template={input.silva_ref})
        screen.seqs(fasta={output.align}, start=1046, end=43116)
        filter.seqs(fasta={output.screen}, trump=., vertical=T)
        unique.seqs(fasta={output.filter})
        pre.cluster(fasta={output.unique}, name={output.names}, diffs=2)
        '
        """

rule calc_dists_rdp_full:
    input:
        fasta=rules.get_full_length_rdp.output.preclust_fasta
    output:
        dist="data/rdp/rdp.bacteria.good.filter.unique.precluster.dist"
    shell:
        """
        {mothur}
        dist.seqs(fasta={input.fasta}, cutoff=0.03)
        '
        """

rule cluster_rdp_full:
    input:
        dist=rules.calc_dists_rdp_full.output.dist,
        names=rules.get_full_length_rdp.output.preclust_names
    output:
        list="data/rdp/rdp.bacteria.good.filter.unique.precluster.opti_mcc.list",
        steps="data/rdp/rdp.bacteria.good.filter.unique.precluster.opti_mcc.steps",
        sensspec="data/rdp/rdp.bacteria.good.filter.unique.precluster.opti_mcc.sensspec"
    shell:
        """
        {mothur}
        cluster(column={input.dist}, name={input.names}, cutoff=0.03)
        '
        """

rule otu_reps_rdp_full:
    input:
        tax=rules.get_rdp_bact.output.tax,
        list=rules.cluster_rdp_full.output.list,
        names=rules.get_full_length_rdp.output.preclust_names
    output:
        tax="data/rdp/rdp.bacteria.good.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy",
        sum="data/rdp/rdp.bacteria.good.filter.unique.precluster.opti_mcc.0.03.cons.tax.summary",
        names="data/rdp/rdp.bacteria.good.filter.unique.precluster.opti_mcc.0.03.rep.names"
    shell:
        """
        {mothur}
        classify.otu(taxonomy={input.tax}, list={input.list}, name={input.names})
        get.oturep(method=abundance, list={input.list}, name={input.names})
        '
        """

rule get_rdp_v4:
    input:
        fasta=rules.get_rdp_bact.output.fasta
    output:
        fasta="data/rdp/rdp.bact.v4.fasta",
        filter="data/rdp/rdp.bact.v4.filter.fasta",
        names="data/rdp/rdp.bact.v4.filter.names",
        unique="data/rdp/rdp.bact.v4.filter.unique.fasta",
        preclust_fasta="data/rdp/rdp.bact.v4.filter.unique.precluster.fasta",
        preclust_names="data/rdp/rdp.bact.v4.filter.unique.precluster.names"
    shell:
        """
        {mothur}
        pcr.seqs(fasta={input.fasta}, start=11894, end=25319, keepdots=F) '
        mv data/rdp/rdp.bacteria.pcr.fasta {output.fasta}
        {mothur}
        filter.seqs(fasta={output.fasta}, trump=., vertical=T)
        unique.seqs(fasta={output.filter})
        pre.cluster(fasta={output.unique}, name={output.names}, diffs=2)
        '
        """

rule calc_dists_rdp_v4:
    input:
        fasta=rules.get_rdp_v4.output.preclust_fasta
    output:
        dist="data/rdp/rdp.bact.v4.filter.unique.precluster.dist"
    shell:
        """
        {mothur}
        dist.seqs(fasta={input.fasta}, cutoff=0.03)
        '
        """

rule cluster_rdp_v4:
    input:
        dist=rules.calc_dists_rdp_v4.output.dist,
        names=rules.get_rdp_v4.output.preclust_names
    output:
        list="data/rdp/rdp.bact.v4.filter.unique.precluster.opti_mcc.list",
        steps="data/rdp/rdp.bact.v4.filter.unique.precluster.opti_mcc.steps",
        sensspec="data/rdp/rdp.bact.v4.filter.unique.precluster.opti_mcc.sensspec"
    shell:
        """
        {mothur}
        cluster(column={input.dist}, name={input.names}, cutoff=0.03)
        '
        """

rule otu_reps_rdp_v4:
    input:
        tax=rules.get_rdp_bact.output.tax,
        list=rules.cluster_rdp_v4.output.list,
        names=rules.get_rdp_v4.output.preclust_names
    output:
        tax="data/rdp/rdp.bact.v4.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy",
        sum="data/rdp/rdp.bact.v4.filter.unique.precluster.opti_mcc.0.03.cons.tax.summary",
        names="data/rdp/rdp.bact.v4.filter.unique.precluster.opti_mcc.0.03.rep.names"
    shell:
        """
        {mothur}
        classify.otu(taxonomy={input.tax}, list={input.list}, name={input.names})
        get.oturep(method=abundance, list={input.list}, name={input.names})
        '
        """
