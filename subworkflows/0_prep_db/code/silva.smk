" Download the reference databases and process with mothur "
import os
import subprocess

wildcard_constraints:
    workdir="data/silva/"

rule silva_targets:
    input:
        "data/silva/silva.bacteria.good.filter.fasta",
        "data/silva/silva.v4.tax"

rule download_silva_db:
    output:
        "{workdir}Silva.seed_v132.tgz"
    shell:
        'wget -N -P {workdir} http://www.mothur.org/w/images/3/32/Silva.seed_v132.tgz'

rule unpack_silva_db:
    input:
        tar="{workdir}Silva.seed_v132.tgz"
    output:
        "{workdir}silva.seed_v132.align",
        "{workdir}silva.seed_v132.tax"
    shell:
        "tar xvzf {input.tar} -C {workdir}"

rule get_silva_prok_lineage:
    input:
        fasta="{workdir}silva.seed_v132.align",
        tax="{workdir}silva.seed_v132.tax"
    output:
        fasta="{workdir}silva.bact_archaea.fasta",
        tax="{workdir}silva.bact_archaea.tax"
    params:
        mothur=mothur_bin
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria-Archaea)' ; "
        "mv {workdir}silva.nr_v132.pick.align {output.fasta} ; "
        "mv {workdir}silva.nr_v132.pick.tax {output.tax}"

rule get_silva_bact_lineage:
    input:
        fasta='{workdir}silva.bact_archaea.fasta',
        tax='{workdir}silva.bact_archaea.tax'
    output:
        fasta="{workdir}silva.bacteria.fasta"
        tax="{workdir}silva.bacteria.tax",
        sum="{workdir}silva.bacteria.summary"
    params:
        mothur=mothur_bin,
        pick_fasta='{workdir}silva.bact_archaea.pick.fasta',
        pick_tax='{workdir}silva.bact_archaea.pick.tax'
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria)'; "
        "mv {params.pick_fasta} {output.fasta}; "
        "mv {params.pick_tax} {output.tax}; "
        "{params.mothur} '#summary.seqs(fasta={output.fasta})'"

rule get_full_length:
    input:
        "{workdir}silva.bacteria.fasta"
    output:
        fasta="{workdir}silva.bacteria.good.fasta",
        sum="{workdir}silva.bacteria.good.summary",
        bad="{workdir}silva.bacteria.bad.accnos"
    params:
        mothur=mothur_bin
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "screen.seqs(fasta={input}, start=1046, end=43116); "
        "summary.seqs(fasta={output.fasta})'"

rule filter_full_length:
    input:
        fasta="{workdir}silva.bacteria.good.fasta"
    output:
        "{workdir}silva.bacteria.good.filter.unique.precluster.fasta",
        "{workdir}silva.bacteria.good.filter.unique.precluster.names",
        "{workdir}silva.bacteria.good.filter.unique.precluster.map"
    params:
        mothur=mothur_bin,
        filter="{workdir}silva.bacteria.good.filter.fasta",
        names="{workdir}silva.bacteria.good.filter.names","
        unique="{workdir}silva.bacteria.good.filter.unique.fasta
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "filter.seqs(fasta={input.fasta}, trump=., vertical=T); "
        "unique.seqs(fasta={params.filter}); "
        "pre.cluster(fasta={params.unique}, name={params.names}, diffs=10)'"

rule trim_v4:
    input:
        fasta="{workdir}silva.bacteria.fasta"
    output:
        "{workdir}silva.bacteria.v4.fastda"
    params:
        mothur=mothur_bin,
        pcr="{workdir}silva.bacteria.pcr.fasta"
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "pcr.seqs(fasta={input.fasta}, start=118994, end=25319, keepdots=F); "
        "summary.seqs(fasta={output}); '"
        "mv {params.pcr} {output}"

rule filter_v4:
    input:
        "{workdir}silva.bacteria.v4.align"
    output:
        "{workdir}silva.bacteria.v4.filter.unique.precluster.fasta",
        "{workdir}silva.bacteria.v4.filter.unique.precluster.names",
        "{workdir}silva.bacteria.v4.filter.unique.precluster.map"
    params:
        mothur=mothur_bin,
        filter="{workdir}silva.bacteria.v4.filter.fasta",
        names="{workdir}silva.bacteria.v4.filter.names","
        unique="{workdir}silva.bacteria.v4.filter.unique.fasta
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "filter.seqs(fasta={input.fasta}, trump=., vertical=T); "
        "unique.seqs(fasta={params.filter}); "
        "pre.cluster(fasta={params.unique}, name={params.names}, diffs=2)'"

rule calc_dists:
    input:
        "{workdir}silva.bacteria.{subset}.filter.unique.precluster.fasta"
    output:
        "{worddir}silva.bacteria.{subset}.filter.unique.precluster.dist",
        "{worddir}silva.bacteria.{subset}.filter.unique.precluster.names"
    wildcard_constraints:
        subset="v4|good"
    params:
        mothur=mothur_bin
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "dist.seqs(fasta={input}, cutoff=0.03)'"

rule cluster:
    input:
        dist="{worddir}silva.bacteria.{subset}.filter.unique.precluster.dist",
        names="{worddir}silva.bacteria.{subset}.filter.unique.precluster.names"
    output:
        list="{worddir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.list",
        steps="{worddir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.steps",
        sensspec="{worddir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.sensspec"
    params:
        mothur=mothur_bin
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "'cluster(column={input.dist}, name={input.names}, cutoff=0.03)'"

rule get_otu_reps:
    input:
        tax="{workdir}silva.bacteria.tax",
        list="{worddir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.list",
        names="{worddir}silva.bacteria.{subset}.filter.unique.precluster.names"
    output:
        "{workdir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy",
        "{workdir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.cons.tax.summary",
        "{workdir}silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.rep.names"
    params:
        mothur=mothur_bin
    shell:
        "{params.mothur} '#set.dir(input={workdir}, output={workdir}); "
        "classify.otu(taxonomy={input.tax}, list={input.list}, names={input.names})"
        "get.oturep(method=abundance, list={input.list}, name={input.names})'"
