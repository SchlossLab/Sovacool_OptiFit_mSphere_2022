
rule remove_gaps_query:
    input:
        fna=prep_samples("data/{dataset}/processed/{dataset}.fasta")
    output:
        fna="data/{dataset}/{dataset}.fna"
    shell:
        """
        cat {input.fna} | sed 's/-//g' > {output.fna}
        """

rule import_fasta:
    input:
        fna=rules.remove_gaps_query.output.fna
    output:
        qza="data/{dataset}/{dataset}_seqs.qza"
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime tools import \
            --input-path {input.fna} \
            --output-path {output.qza} \
            --type 'FeatureData[Sequence]'
        """


rule prep_count_table_biom:
    input:
        count=prep_samples("data/{dataset}/processed/{dataset}.count_table"),
        code="code/prep_count_table.R"
    output:
        count='data/{dataset}/{dataset}.count_table.tsv',
    script:
        "code/prep_count_table.R"



rule make_biom:
    input:
        count=rules.prep_count_table_biom.output.count
    output:
        biom='data/{dataset}/{dataset}.biom'
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        biom convert -i {input.count} -o {output.biom} --table-type="OTU table" --to-hdf5
        """

rule import_table:
    input:
        biom=rules.make_biom.output.biom
    output:
        qza="data/{dataset}/{dataset}_table.qza"
    params:
        outdir='data/{dataset}/'
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime tools import \
            --input-path {input.biom} \
            --type 'FeatureTable[Frequency]' \
            --input-format BIOMV210Format \
            --output-path {output.qza}
        """

rule remove_gaps_ref:
    input:
        fna=fit_ref_db("data/{dataset}_{ref}_{region}/preclust_db/{ref}.{region}.filter.unique.precluster.fasta")
    output:
        fna="data/{dataset}/{ref}_{region}.fna"
    shell:
        """
        cat {input.fna} | sed 's/-//g' | sed 's/^[^>]*\.//g' > {output.fna}
        """

rule import_ref:
    input:
        fna=rules.remove_gaps_ref.output.fna
    output:
        qza="data/{dataset}/{ref}_{region}_seqs.qza"
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime tools import \
            --input-path {input.fna} \
            --output-path {output.qza} \
            --type 'FeatureData[Sequence]'
        """

rule qiime_de_novo:
    input:
        table=rules.import_table.output.qza,
        seqs=rules.import_fasta.output.qza
    output:
        table='results/{dataset}/de_novo/table.qza',
        seqs='results/{dataset}/de_novo/seqs.qza'
    benchmark:
        'benchmarks/{dataset}/qiime_de_novo.txt'
    params:
        p=perc_identity
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime vsearch cluster-features-de-novo \
            --i-table {input.table} \
            --i-sequences {input.seqs} \
            --p-perc-identity {params.p} \
            --o-clustered-table {output.table} \
            --o-clustered-sequences {output.seqs}
        """

rule qiime_closed_ref:
    input:
        table=rules.import_table.output.qza,
        seqs=rules.import_fasta.output.qza,
        ref=rules.import_ref.output.qza
    output:
        table='results/{dataset}/closed_ref/{ref}/{region}/table.qza',
        seqs='results/{dataset}/closed_ref/{ref}/{region}/seqs.qza',
        unmatched='results/{dataset}/closed_ref/{ref}/{region}/unmatched.qza'
    benchmark:
        'benchmarks/{dataset}/closed_ref.{ref}.{region}.txt'
    params:
        p=perc_identity
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime vsearch cluster-features-closed-reference \
            --i-table {input.table} \
            --i-sequences {input.seqs} \
            --i-reference-sequences {input.ref} \
            --p-perc-identity {params.p} \
            --o-clustered-table {output.table} \
            --o-clustered-sequences {output.seqs} \
            --o-unmatched-sequences {output.unmatched}
        """

rule qiime_open_ref:
    input:
        table=rules.import_table.output.qza,
        seqs=rules.import_fasta.output.qza,
        ref=rules.import_ref.output.qza
    output:
        table='results/{dataset}/open_ref/{ref}/{region}/table.qza',
        seqs='results/{dataset}/open_ref/{ref}/{region}/seqs.qza',
        new='results/{dataset}/open_ref/{ref}/{region}/new.qza'
    benchmark:
        'benchmarks/{dataset}/open_ref.{ref}.{region}.txt'
    params:
        p=perc_identity
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime vsearch cluster-features-open-reference \
            --i-table {input.table} \
            --i-sequences {input.seqs} \
            --i-reference-sequences {input.ref} \
            --p-perc-identity {params.p} \
            --o-clustered-table {output.table} \
            --o-clustered-sequences {output.seqs} \
            --o-new-reference-sequences {output.new}
        """

rule export_qiime:
    input:
        table="{input_dir}/table.qza",
        seqs="{input_dir}/seqs.qza"
    output:
        biom="{input_dir}/feature-table.biom",
        fna="{input_dir}/dna-sequences.fasta"
    params:
        expdir="{input_dir}/"
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        qiime tools export \
            --input-path {input.table} \
            --output-path {params.expdir}
        qiime tools export \
            --input-path {input.seqs} \
            --output-path {params.expdir}
        """

rule biom_to_tsv:
    input:
        biom=rules.export_qiime.output.biom
    output:
        tsv="{input_dir}/feature-table.tsv"
    conda:
        'config/qiime2-2020.8-py36-linux-conda.yml'
    shell:
        """
        biom convert -i {input.biom} -o {output.tsv} --to-tsv
        """


rule biom_tsv_to_list:
    input:
        R="code/biom_tsv_to_list.py",
        tsv=rules.biom_to_tsv.output.tsv
    output:
        list="LISTFILE"
    params:
        perc_identity=perc_identity # label = 1 - perc_identity
    script:
        "code/biom_tsv_to_list.py"