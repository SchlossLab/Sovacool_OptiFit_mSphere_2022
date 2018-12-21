""" Benchmarking the OptiFit algorithm using an external reference database """

rule fit_to_external_ref_db:
    input:
        expand('results/{{reference}}-as-reference/{{dataset}}/{{dataset}}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
        expand('results/{{reference}}-as-reference/{{dataset}}/{{dataset}}_i-{iter}/r-{rep}/method-{method}_printref-f/sample.optifit_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods),
        "results/{reference}-as-reference/{dataset}/aggregate.sensspec"

rule prep_reference:
    input:
        accnos="data/references/{reference}.accnos",
        fasta="data/references/{reference}.align",
        column="data/references/{reference}.dist"
    output:
        accnos="results/{reference}-as-reference/{reference}.accnos",
        fasta="results/{reference}-as-reference/{reference}.fasta",
        name="results/{reference}-as-reference/{reference}.names",
        unique="results/{reference}-as-reference/{reference}.unique.fasta"
    params:
        output_dir="{output_dir}/{reference}-as-reference/"
    benchmark:
        "benchmarks/{reference}-as-reference/{reference}.prep_reference.log"
    log:
        "logfiles/{reference}-as-reference/{reference}.prep_reference.log"
    shell:
        'cp {input.accnos} {output.accnos}; cp {input.fasta} {output.fasta}; '
        '{params.mothur} "#set.logfile(name={log}); set.dir(output={params.output_dir}); unique.seqs(fasta={input.fasta})" '

rule reference_cluster:
    input:
        names="data/references/{reference}/{reference}.names",
        column="data/references/{reference}/{reference}.dist"
    output:
        expand('results/{{reference}}-as-reference/{{reference}}.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
    params:
        mothur=mothur_bin,
        output_dir='results/{reference}-as-reference/'
    benchmark:
        "benchmarks/{reference}-as-reference/{reference}.cluster.log"
    log:
        "logfiles/{reference}-as-reference/{reference}.cluster.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.rep}); set.dir(output={params.output_dir}); cluster(column={input.column}, name={input.names}, cutoff=0.3)"'


rule sample_cluster:
    input:
        count="data/{dataset}/{dataset}.count_table",
        column="data/{dataset}/{dataset}.dist"
    output:
        expand('results/{{reference}}-as-reference/{{dataset}}/{{dataset}}__i-{{iter}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
    params:
        mothur=mothur_bin,
        output_dir='results/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/',
        iter="{iter}"
    benchmark:
        "benchmarks/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/sample.cluster.log"
    log:
        "logfiles/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/sample.cluster.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.iter}); set.dir(output={params.output_dir}); cluster(column={input.column}, count={input.count}, cutoff=0.3)"'

rule fit_to_ref:
    input:
        reflist='results/{reference}-as-reference/{reference}.opti_mcc.list',
        refcolumn='results/{reference}-as-reference/{reference}.dist',
        refcount='results/{reference}-as-reference/{reference}.count_table',
        reffasta='results/{reference}-as-reference/{reference}.fasta',
        fasta='results/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/sample/sample.fasta',
        count='results/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/sample/sample.count_table',
        column='results/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/sample/sample.dist'
    output:
        expand('results/{{reference}}-as-reference/{{dataset}}/{{dataset}}__i-{{iter}}/r-{{rep}}/method-{{method}}_printref-f/sample.optifit_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
    params:
        mothur=mothur_bin,
        output_dir="results/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/r-{rep}/method-{method}_printref-f/",
        rep="{rep}",
        method="{method}"
    benchmark:
        "benchmarks/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/r-{rep}/method-{method}_printref-f/fit.log"
    log:
        "logfiles/{reference}-as-reference/{dataset}/{dataset}_i-{iter}/r-{rep}/method-{method}_printref-f/fit.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.rep}); set.dir(output={params.output_dir}); cluster.fit(reflist={input.reflist}, refcolumn={input.refcolumn}, refcount={input.refcount}, reffasta={input.reffasta}, fasta={input.fasta}, count={input.count}, column={input.column}, printref=f, method={params.method})"'

rule ref_aggregate_sensspec:
    input:
        opticlust=expand('results/{{reference}}-as-reference/{{dataset}}/{{dataset}}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
        optifit=expand('results/{{reference}}-as-reference/{{dataset}}/{{dataset}}_i-{iter}/r-{rep}/method-{method}_printref-f/sample.optifit_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods)
    output:
        "results/{reference}-as-reference/{dataset}/aggregate.sensspec"
    params:
        dataset="{dataset}",
        reference_fractions=reference_fractions,
        weights=weights,
        iters=iters,
        reps=reps,
        methods=methods,
        prefixes=['sample','reference']
    benchmark:
        "benchmarks/{reference}-as-reference/{dataset}/aggregate_sensspec.log"
    run:
        header_str = 'iter\tlabel\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\treference_fraction\titer\trep\ttype\n'
        with open(output[0], 'w') as output_file:
            output_file.write(header_str)
            for weight in params.weights:
                for reference_fraction in params.reference_fractions:
                    for iter in params.iters:
                        for rep in params.reps:
                            for prefix in params.prefixes:
                                input_filename = f'results/{reference}-as-reference/{params.dataset}/{params.dataset}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec'
                                with open(input_filename, 'r') as input_file:
                                    for line in input_file:
                                        pass
                                    opticlust_result = re.sub(r'(\S*\t\S*\t)(.*)', r'\t\1\t\2', line).rstrip()
                                    output_file.write(f"{opticlust_result}\t{reference_fraction}\t{iter}\t{rep}\t{prefix}\n")
                            for method in params.methods:
                                input_filename = f"results/{reference}-as-reference/{params.dataset}/{params.dataset}_i-{iter}/r-{rep}/method-{method}_printref-f/sample.optifit_mcc.sensspec"
                                with open(input_filename, 'r') as input_file:
                                    for line in input_file:
                                        pass
                                    line = line.strip()
                                    output_file.write(f"{line}\t{reference_fraction}\t{iter}\t{rep}\tmethod-{method}_printref-f\n")
