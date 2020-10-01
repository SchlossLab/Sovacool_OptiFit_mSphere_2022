
rule merge_sensspec1:
    priority: 5
    input:
        code='../../code/R/merge_results.R',
        fcns='../../code/R/functions.R',
        tsv=expand('results/{dataset}/refweight_{ref_weight}/fitpercent_{fitpercent}/seed_{seed}/method_{method}/printref_{printref}/{dataset}.pick.optifit_mcc.mod.sensspec',
            dataset=datasets, ref_weight=random_refweight_options, fitpercent=fitpercents, seed=seeds, method=methods, printref=printrefs)
    output:
        tsv="results/pre-join/sensspec.tsv"
    log:
        'log/merge_sensspec.txt'
    benchmark:
        'benchmarks/merge_sensspec.txt'
    script:
        '../../../code/R/merge_results.R'

rule merge_benchmarks1:
    priority: 7
    input:
        code='../../code/R/merge_results.R',
        fcns='../../code/R/functions.R',
        tsv=expand('benchmarks/{dataset}/fit.method_{method}.printref_{printref}.refweight_{ref_weight}.fitpercent_{fitpercent}.seed_{seed}.mod.txt',
            dataset=datasets, method=methods, printref=printrefs, ref_weight=random_refweight_options, fitpercent=fitpercents, seed=seeds)
    output:
        tsv="results/pre-join/benchmarks.tsv"
    log:
        'log/merge_benchmarks.txt'
    benchmark:
        'benchmarks/merge_benchmarks.txt'
    script:
        '../../../code/R/merge_results.R'

rule merge_fractions1:
    priority: 1
    input:
        code='../../code/R/merge_results.R',
        fcns='../../code/R/functions.R',
        tsv=expand('results/{dataset}/refweight_{ref_weight}/fitpercent_{fitpercent}/seed_{seed}/method_{method}/printref_{printref}/fraction_reads_mapped.txt',
            dataset=datasets, method=['closed'], printref=printrefs, ref_weight=random_refweight_options, fitpercent=fitpercents, seed=seeds)
    output:
        tsv="results/pre-join/fraction_reads_mapped.tsv"
    log:
        'log/merge_fractions.txt'
    benchmark:
        'benchmarks/merge_fractions.txt'
    script:
        '../../../code/R/merge_results.R'

rule merge_sizes1:
    input:
        code='../../code/R/merge_results.R',
        fcns='../../code/R/functions.R',
        tsv=expand("results/{dataset}/refweight_{ref_weight}/fitpercent_{fitpercent}/seed_{seed}/input_size.tsv",
            dataset=datasets, ref_weight=random_refweight_options, fitpercent=fitpercents, seed=seeds)
    output:
        tsv="results/pre-join/input_sizes.tsv"
    log:
        'log/merge_sizes.txt'
    benchmark:
        'benchmarks/merge_sizes.txt'
    script:
        '../../../code/R/merge_results.R'

rule merge_diversity1:
    priority: 3
    input:
        code='../../code/R/merge_results.R',
        fcns='../../code/R/functions.R',
        tsv=expand('results/{dataset}/refweight_{ref_weight}/fitpercent_{fitpercent}/seed_{seed}/method_{method}/printref_{printref}/{dataset}.pick.optifit_mcc.summary',
            dataset=datasets, ref_weight=random_refweight_options, fitpercent=fitpercents, seed=seeds, method=methods, printref=printrefs)
    output:
        tsv='results/pre-join/diversity.tsv'
    log:
        'log/merge_gap_counts.txt'
    script:
        '../../../code/R/merge_results.R'
