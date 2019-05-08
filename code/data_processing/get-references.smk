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
        f"{output_dir}/silva/Silva.seed_{version}.tgz"
    shell:
        f'wget -N -P {output_dir}/silva/ http://www.mothur.org/w/images/3/32/Silva.seed_{version}.tgz'

rule unpack_silva_db:
    input:
        tar=f"{output_dir}/silva/Silva.seed_{version}.tgz"
    output:
        f"{output_dir}/silva/silva.seed_{version}.align",
        f"{output_dir}/silva/silva.seed_{version}.tax"
    shell:
        f"tar xvzf {{input.tar}} -C {output_dir}/"

rule get_prok_lineage:
    input:
        fasta=f"{output_dir}/silva/silva.seed_{version}.align",
        tax=f"{output_dir}/silva/silva.seed_{version}.tax"
    output:
        fasta=f"{output_dir}/silva/silva.bact_archaea.align",
        tax=f"{output_dir}/silva/silva.bact_archaea.tax"
    params:
        mothur=mothur_bin,
        version=version
    shell:
        "{params.mothur} '#get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria-Archaea)' ; "
        "mv {output_dir}/silva/silva.nr_{params.version}.pick.align {output.fasta} ; "
        "mv {output_dir}/silva/silva.nr_{params.version}.pick.tax {output.tax}"

rule get_bact_lineage:
    input:
        fasta=f'{output_dir}/silva/silva.bact_archaea.align',
        tax=f'{output_dir}/silva/silva.bact_archaea.tax'
    output:
        f'{output_dir}/silva/silva.bact_archaea.pick.align',
        f'{output_dir}/silva/silva.bact_archaea.pick.tax'
    params:
        mothur=mothur_bin,
    shell:
        "{params.mothur} '#get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria)'"

rule rename_bact:
    input:
        "{output_dir}/silva/silva.bact_archaea.pick.{ext}"
    output:
        "{output_dir}/silva/silva.bacteria.{ext}"
    wildcard_constraints:
        ext="align|tax"
    shell:
        "mv {input} {output}"

rule pcr_seqs:
    input:
        "{output_dir}/silva/silva.bacteria.align"
    output:
        "{output_dir}/silva/silva.bacteria.pcr.ng.names",
        "{output_dir}/silva/silva.bacteria.pcr.align"
    params:
        mothur=mothur_bin,
    shell:
        '{params.mothur} "#pcr.seqs(fasta={input}, start=13862, end=23445, keepdots=F);degap.seqs();unique.seqs()"'

rule get_pcr_accession_numbers:
    input:
        "{output_dir}/silva/silva.bacteria.pcr.ng.names"
    output:
        "{output_dir}/silva/silva.bacteria.pcr.ng.accnos"
    shell:
        'cut -f 1 {input} > {output}'

rule get_fasta_seqs:
    input:
        fasta="{output_dir}/silva/silva.bacteria.pcr.align",
        accnos="{output_dir}/silva/silva.bacteria.pcr.ng.accnos"
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
        "{output_dir}/silva/silva.v4.align"
    output:
        "{output_dir}/silva/silva.v4.accnos"
    shell:
        'grep "^>" {input} | cut -c 2- > {output}'

rule get_taxon_seqs:
    input:
        tax="{output_dir}/silva/silva.bacteria.tax",
        accnos="{output_dir}/silva/silva.v4.accnos"
    output:
        tax="{output_dir}/silva/silva.v4.tax"
    params:
        mothur=mothur_bin,
        pick="{output_dir}/silva/silva.bacteria.pick.tax"
    shell:
        '{params.mothur} "#get.seqs(taxonomy={input.tax}, accnos={input.accnos})"; '
        'mv {params.pick} {output.tax}'

rule download_rdp_db:
    output:
        expand("{{output_dir}}/rdp/trainset14_032015.pds/trainset14_032015.pds.{ext}",ext={'tax', 'fasta'})
    params:
        dir="{output_dir}/rdp/",
        tar="Trainset14_032015.pds.tgz"
    shell:
        "wget -N -P {params.dir} http://www.mothur.org/w/images/8/88/{params.tar} ; "
        "tar xvzf {params.dir}{params.tar} -C {params.dir} ; "
