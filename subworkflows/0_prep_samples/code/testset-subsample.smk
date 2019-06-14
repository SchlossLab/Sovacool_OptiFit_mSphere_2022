import os

samples = config['samples']
subsample_size = config['subsample_size']

rule subset_all:
    input:
        expand("data/processed/{sample}_{subsample_size}/{sample}_{subsample_size}.{ext}", input_dir=input_dir, subsample_size=subsample_size, sample=samples, ext={'fasta', 'accnos', 'count_table'})

rule subset:
    input:
        f'data/raw/{{sample}}.fasta'
    output:
        expand("data/raw/{{sample}}_{{subsample_size}}/{{sample}}_{{subsample_size}}.{ext}", input_dir=input_dir, ext={'fasta', 'accnos', 'count_table'})
    params:
        mothur=mothur_bin,
        sample = '{sample}',
        input_dir = 'data/raw/',
        output_dir = 'data/processed/{{sample}}_{{subsample_size}}/'
    shell:
        '{params.mothur} "#set.dir(input={params.input_dir}{params.sample}, output={params.output_dir}); '
        'sub.sample(fasta={params.sample}.fasta, size={subsample_size}); '
        'list.seqs(fasta=current); '
        'get.seqs(accnos=current, count={params.sample}.count_table); '
        'rename.file(accnos = current, fasta=current, count=current, prefix={params.sample}_{subsample_size});"'
