" Download the reference databases and process with mothur "
import os
import subprocess

rule silva_targets:
    input:
        expand("data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.{ext}", subset=("v4", "good"), ext=("cons.taxonomy", "cons.tax.summary", "rep.names"))

rule download_silva_db:
    output:
        temp("data/silva/Silva.seed_v132.tgz")
    shell:
        'wget -N -P data/silva/ http://www.mothur.org/w/images/3/32/Silva.seed_v132.tgz'

rule unpack_silva_db:
    input:
        tar="data/silva/Silva.seed_v132.tgz"
    output:
        "data/silva/silva.seed_v132.align",
        "data/silva/silva.seed_v132.tax"
    shell:
        "tar xvzf {input.tar} -C data/silva/"

rule get_silva_prok_lineage:
    input:
        fasta="data/silva/silva.seed_v132.align",
        tax="data/silva/silva.seed_v132.tax"
    output:
        fasta="data/silva/silva.bact_archaea.fasta",
        tax="data/silva/silva.bact_archaea.tax"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria-Archaea)' ; "
        "mv data/silva/silva.nr_v132.pick.align {output.fasta} ; "
        "mv data/silva/silva.nr_v132.pick.tax {output.tax}"

rule get_silva_bact_lineage:
    input:
        fasta='data/silva/silva.bact_archaea.fasta',
        tax='data/silva/silva.bact_archaea.tax'
    output:
        fasta="data/silva/silva.bacteria.fasta",
        tax="data/silva/silva.bacteria.tax",
        sum="data/silva/silva.bacteria.summary"
    params:
        pick_fasta='data/silva/silva.bact_archaea.pick.fasta',
        pick_tax='data/silva/silva.bact_archaea.pick.tax'
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria)'; "
        "mv {params.pick_fasta} {output.fasta}; "
        "mv {params.pick_tax} {output.tax}; "
        "mothur '#summary.seqs(fasta={output.fasta})'"

rule get_full_length:
    input:
        "data/silva/silva.bacteria.fasta"
    output:
        fasta="data/silva/silva.bacteria.good.fasta",
        sum="data/silva/silva.bacteria.good.summary",
        bad="data/silva/silva.bacteria.bad.accnos"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "screen.seqs(fasta={input}, start=1046, end=43116); "
        "summary.seqs(fasta={output.fasta})'"

rule filter_full_length:
    input:
        fasta="data/silva/silva.bacteria.good.fasta"
    output:
        "data/silva/silva.bacteria.good.filter.unique.precluster.fasta",
        "data/silva/silva.bacteria.good.filter.unique.precluster.names",
        "data/silva/silva.bacteria.good.filter.unique.precluster.map"
    params:
        filter="data/silva/silva.bacteria.good.filter.fasta",
        names="data/silva/silva.bacteria.good.filter.names",
        unique="data/silva/silva.bacteria.good.filter.unique.fasta"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "filter.seqs(fasta={input.fasta}, trump=., vertical=T); "
        "unique.seqs(fasta={params.filter}); "
        "pre.cluster(fasta={params.unique}, name={params.names}, diffs=10)'"

rule trim_v4:
    input:
        fasta="data/silva/silva.bacteria.fasta"
    output:
        "data/silva/silva.bacteria.v4.fasta"
    params:
        pcr="data/silva/silva.bacteria.pcr.fasta"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "pcr.seqs(fasta={input.fasta}, start=118994, end=25319, keepdots=F); "
        "summary.seqs(fasta={output}); '"
        "mv {params.pcr} {output}"

rule filter_v4:
    input:
        fasta="data/silva/silva.bacteria.v4.fasta"
    output:
        "data/silva/silva.bacteria.v4.filter.unique.precluster.fasta",
        "data/silva/silva.bacteria.v4.filter.unique.precluster.names",
        "data/silva/silva.bacteria.v4.filter.unique.precluster.map"
    params:
        filter="data/silva/silva.bacteria.v4.filter.fasta",
        names="data/silva/silva.bacteria.v4.filter.names",
        unique="data/silva/silva.bacteria.v4.filter.unique.fasta"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "filter.seqs(fasta={input.fasta}, trump=., vertical=T); "
        "unique.seqs(fasta={params.filter}); "
        "pre.cluster(fasta={params.unique}, name={params.names}, diffs=2)'"

rule calc_dists:
    input:
        "data/silva/silva.bacteria.{subset}.filter.unique.precluster.fasta"
    output:
        "data/silva/silva.bacteria.{subset}.filter.unique.precluster.dist"
    wildcard_constraints:
        subset="v4|good"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "dist.seqs(fasta={input}, cutoff=0.03)'"

rule cluster:
    input:
        dist="data/silva/silva.bacteria.{subset}.filter.unique.precluster.dist",
        names="data/silva/silva.bacteria.{subset}.filter.unique.precluster.names"
    output:
        list="data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.list",
        steps="data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.steps",
        sensspec="data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.sensspec"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "'cluster(column={input.dist}, name={input.names}, cutoff=0.03)'"

rule get_otu_reps:
    input:
        tax="data/silva/silva.bacteria.tax",
        list="data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.list",
        names="data/silva/silva.bacteria.{subset}.filter.unique.precluster.names"
    output:
        "data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy",
        "data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.cons.tax.summary",
        "data/silva/silva.bacteria.{subset}.filter.unique.precluster.opti_mcc.0.03.rep.names"
    shell:
        "mothur '#set.dir(input=data/silva/, output=data/silva/); "
        "classify.otu(taxonomy={input.tax}, list={input.list}, names={input.names}); "
        "get.oturep(method=abundance, list={input.list}, name={input.names})'"
