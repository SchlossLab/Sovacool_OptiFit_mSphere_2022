" Download the reference database and process with mothur "
import os
import subprocess

# TODO: add greengenes download

output_dir = os.path.join(config['input_dir'], 'references')
version = config['silva_db_version']

rule ref_db_targets:
    input:
        f"{output_dir}/silva/silva.v4.tax",
        expand("{output_dir}/rdp/trainset14_032015.pds/trainset14_032015.pds.{ext}", output_dir=output_dir, ext={'tax', 'fasta'})

rule download_silva_db:
    output:
        f"data/references/silva/Silva.seed_{version}.tgz"
    shell:
        f'wget -N -P data/references/silva/ http://www.mothur.org/w/images/3/32/Silva.seed_{version}.tgz'

rule unpack_silva_db:
    input:
        tar=f"data/references/silva/Silva.seed_{version}.tgz"
    output:
        f"data/references/silva/silva.seed_{version}.align",
        f"data/references/silva/silva.seed_{version}.tax"
    shell:
        "tar xvzf {input.tar} -C data/references/"

rule get_prok_lineage:
    input:
        fasta=f"data/references/silva/silva.seed_{version}.align",
        tax=f"data/references/silva/silva.seed_{version}.tax"
    output:
        fasta="data/references/silva/silva.bact_archaea.fasta",
        tax="data/references/silva/silva.bact_archaea.tax"
    params:
        mothur=mothur_bin,
        version=version
    shell:
        "{params.mothur} '#get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria-Archaea)' ; "
        "mv data/references/silva/silva.nr_{params.version}.pick.align {output.fasta} ; "
        "mv data/references/silva/silva.nr_{params.version}.pick.tax {output.tax}"

rule get_bact_lineage:
    input:
        fasta='data/references/silva/silva.bact_archaea.fasta',
        tax='data/references/silva/silva.bact_archaea.tax'
    output:
        'data/references/silva/silva.bact_archaea.pick.fasta',
        'data/references/silva/silva.bact_archaea.pick.tax'
    params:
        mothur=mothur_bin,
    shell:
        "{params.mothur} '#get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria)'"

rule rename_bact:
    input:
        "data/references/silva/silva.bact_archaea.pick.{ext}"
    output:
        "data/references/silva/silva.bacteria.{ext}"
    wildcard_constraints:
        ext="fasta|tax"
    shell:
        "mv {input} {output}"

rule pcr_seqs:
    input:
        "data/references/silva/silva.bacteria.fasta"
    output:
        "data/references/silva/silva.bacteria.pcr.ng.names",
        "data/references/silva/silva.bacteria.pcr.fasta"
    params:
        mothur=mothur_bin,
    shell:
        '{params.mothur} "#pcr.seqs(fasta={input}, start=13862, end=23445, keepdots=F);degap.seqs();unique.seqs()"'

rule get_pcr_accession_numbers:
    input:
        "data/references/silva/silva.bacteria.pcr.ng.names"
    output:
        "data/references/silva/silva.bacteria.pcr.ng.accnos"
    shell:
        'cut -f 1 {input} > {output}'

rule get_fasta_seqs:
    input:
        fasta="data/references/silva/silva.bacteria.pcr.fasta",
        accnos="data/references/silva/silva.bacteria.pcr.ng.accnos"
    output:
        "{output_dir}/silva/silva.v4.align"
    params:
        mothur=mothur_bin,
        pick="{output_dir}/silva/silva.bacteria.pcr.pick.good.filter.fasta"
    shell:
        '{params.mothur} "#get.seqs(fasta={input.fasta}, accnos={input.accnos});screen.seqs(minlength=240, maxlength=275, maxambig=0, maxhomop=8); filter.seqs(vertical=T)"; '
        'mv {params.pick} {output}'

rule get_filtered_accession_numbers:
    input:
        "data/references/silva/silva.fasta"
    output:
        "data/references/silva/silva.accnos"
    shell:
        'grep "^>" {input} | cut -c 2- > {output}'

rule get_taxon_seqs:
    input:
        tax="data/references/silva/silva.bacteria.tax",
        accnos="data/references/silva/silva.accnos"
    output:
        tax="data/references/silva/silva.tax"
    params:
        mothur=mothur_bin,
        pick="data/references/silva/silva.bacteria.pick.tax"
    shell:
        '{params.mothur} "#get.seqs(taxonomy={input.tax}, accnos={input.accnos})"; '
        'mv {params.pick} {output.tax}'

rule download_rdp_db:
    output:
        expand("data/references/rdp/trainset14_032015.pds/trainset14_032015.pds.{ext}",ext={'tax', 'fasta'})
    params:
        dir="data/references/rdp/",
        tar="Trainset14_032015.pds.tgz"
    shell:
        "wget -N -P {params.dir} http://www.mothur.org/w/images/8/88/{params.tar} ; "
        "tar xvzf {params.dir}{params.tar} -C {params.dir} ; "
