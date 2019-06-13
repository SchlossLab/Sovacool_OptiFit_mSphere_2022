import requests

configfile: "config.yaml"

class Sample:
    def __init__(self, name, sra_list, url):
        self.name = name
        self.sra_list = sra_list
        self.url = url

    def download(self, filename, sra = None):
        url = f"{self.url}{sra}/{sra}.sra" if sra else self.url
        request = requests.get(url)
        with open(filename, 'wb') as file:
            file.write(request.content)

samples = {"human": Sample("human", [f"SRR{num}" for num in range(2143516, 2144136+1)],
						   "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR062/SRR062005/"),
           "marine": Sample("marine", [f"SRR{num}" for num in range(3085688, 701+1)],
							"ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR308/SRR3085$ACC/SRR3085$ACC.sra"),
           "mouse": Sample("mouse", list(),
						   "http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar"),
           "soil": Sample("soil", [f"ERR{num}" for num in range(1039459, 1039477+1)],
						  "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR103/ERR10394${FILE}/ERR10394${FILE}.sra")
           }

rule all:
	input:
		expand("data/raw/human/{sra}_{num}.fastq.gz", sra=config["human_sra"], num=('1', '2')),
		expand("data/raw/human/human.{prefix}.contigs.{suffix}", prefix=("trim","scrap"), suffix=("fasta", "qual"))

rule download_sra:
	output:
		sra="data/raw/{sample}/SRA/{sra}.sra"
	run:
        samples[wildcards.sample].download(output.sra, wildcards.sra)

rule download_mouse:
    output:
        temp(tar="data/raw/mouse/mouse.tar")
    run:
        samples["mouse"].download(output.tar)

rule extract_mouse_fastq:
    input:
        tar="data/raw/mouse/mouse.tar"
    output:
        fastqs="data/raw/mouse/fastq/{id}.fastq.gz"
    shell:
        "tar xvf {input.tar} -C data/raw/mouse/fastq/"

rule extract_fastq:
	input:
		expand("data/raw/human/SRA/{sra}.sra", sra=config["sras"])
	output:
		expand("data/raw/human/fastq/{sra}_1.fastq", sra=config["sras"]),
		expand("data/raw/human/fastq/{sra}_2.fastq", sra=config["sras"])
	shell:
		"fastq-dump --split-files {input} --outdir data/"

rule names_file:
	input:
		expand("data/{sample}_{num}.fastq", sample=config["{dataset}"], num=('1', '2'))
	output:
		"data/human.files"
	script:
		"code/human_names.R"

rule make_contigs:
	input:
		expand("data/{{dataset}}/{sample}_{num}.fastq.gz", sample=[sample for dataset in config for sample in config[dataset]], num=('1', '2')),
		files="data/{dataset}/human.files"
	output:
		expand("data/{{dataset}}/human.{prefix}.contigs.{suffix}", prefix=("trim","scrap"), suffix=("fasta", "qual")),
		expand("data/{{dataset}}/human.contigs.{suffix2}", suffix2=("report", "groups"))
	params:
		dataset = "{dataset}"
    resources:
        procs=8
	shell:
		"mothur '#set.dir(input=data/{params.dataset}/, output=data/{params.dataset}/); make.contigs(file={input.files}, processors={resources.procs})'"

rule screen_seqs_contigs:
	input:
		fasta="data/{dataset}/human.trim.contigs.fasta",
		groups="data/{dataset}/human.contigs.groups"
	output:
		"data/{dataset}/human.trim.contigs.good.fasta",
		"data/{dataset}/human.trim.contigs.bad.sranos",
		"data/{dataset}/human.contigs.good.groups"
	shell:
		"mothur '#set.dir(input=data, output=data); screen.seqs(fasta={input.fasta}, group={input.groups}, maxambig=0, maxlength=275, maxhomop=8);'"

rule unique_seqs_contigs:
	input:
		fasta="data/human.trim.contigs.good.fasta"
	output:
		"data/human.trim.contigs.good.names",
		"data/human.trim.contigs.good.unique.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); unique.seqs(fasta={input.fasta})'"

rule count_seqs:
	input:
		names="data/human.trim.contigs.good.names",
		groups="data/human.contigs.good.groups"
	output:
		"data/human.trim.contigs.good.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); count.seqs(name={input.names}, group={input.groups});'"

rule align_seqs:
	input:
		fasta="data/human.trim.contigs.good.fasta",
		reference="data/references/silva.v4.align"
	output:
		"data/human.trim.contigs.good.unique.align",
		"data/human.trim.contigs.good.unique.align.report",
		"data/human.trim.contigs.good.unique.flip.sranos"
	shell:
		"mothur '#set.dir(input=data, output=data); align.seqs(fasta={input.fasta}, reference={input.reference}, processors=2);'"

rule screen_seqs_silva:
	input:
		count="data/human.trim.contigs.good.count_table",
		fasta="data/human.trim.contigs.good.unique.align"
	output:
		"data/human.trim.contigs.good.unique.good.align",
		"data/human.trim.contigs.good.unique.bad.sranos",
		"data/human.trim.contigs.good.good.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); screen.seqs(fasta={input.fasta}, count={input.count}, start=5, end=860);'"

rule filter_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.good.align"
	output:
		"data/human.filter",
		"data/human.trim.contigs.good.unique.good.filter.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); filter.seqs(fasta={input.fasta}, vertical=T, trump=.);'"
rule unique_seqs_silva:
	input:
		count="data/human.trim.contigs.good.good.count_table",
		fasta="data/human.trim.contigs.good.unique.good.filter.fasta"
	output:
		"data/human.trim.contigs.good.unique.good.filter.count_table",
		"data/human.trim.contigs.good.unique.good.filter.unique.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); unique.seqs(fasta={input.fasta}, count={input.count});'"

rule precluster:
	input:
		count="data/human.trim.contigs.good.unique.good.filter.count_table",
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.fasta"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); pre.cluster(fasta={input.fasta}, count={input.count}, diffs=2);'"

rule chimera_uchime:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.count_table"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.sranos"
	shell:
		"mothur '#set.dir(input=data, output=data); chimera.uchime(fasta={input.fasta}, count={input.count}, dereplicate=T);'"

rule remove_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.fasta",
		sranos="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.sranos"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta"
	shell:
		"mothur '#set.dir(input=data, output=data); remove.seqs(fasta={input.fasta}, sranos={input.sranos});'"

rule classify_seqs:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		trainset_fasta="data/references/trainset14_032015.pds.fasta",
		trainset_tax="data/references/trainset14_032015.pds.tax"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
	shell:
		"mothur '#set.dir(input=data, output=data); classify.seqs(fasta={input.fasta}, count={input.count}, reference={input.trainset_fasta}, taxonomy={input.trainset_tax}, cutoff=80);'"

rule remove_lineage:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table",
		taxonomy="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy"
	output:
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta",
		"data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table"
	shell:
		"mothur '#set.dir(input=data, output=data); remove.lineage(fasta={input.fasta}, count={input.count}, taxonomy={input.taxonomy}, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);'"

rule rename_final_files:
	input:
		fasta="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta",
		count="data/human.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table",
		taxonomy="data/human.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy"
	output:
		"data/human.fasta",
		"data/human.count_table",
		"data/human.taxonomy"
	shell:
		"mv {input.fasta} data/human.fasta; "
		"mv {input.count} data/human.count_table; "
		"mv {input.taxonomy} data/human.taxonomy"

rule dists:
	input:
		fasta="data/human.fasta"
	output:
		"data/human.dist"
	shell:
		"mothur '#set.dir(input=data, output=data); dist.seqs(fasta={input.fasta}, cutoff=0.03)'"

rule mothur_process:
    shell:
        'mothur "#set.dir(output=$HUMAN);
        	make.contigs(inputdir=$HUMAN, file=human.files, processors=12);
        	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8);
        	unique.seqs();
        	count.seqs(name=current, group=current);
        	align.seqs(fasta=current, reference=data/references/silva.v4.align, processors=2);
        	screen.seqs(fasta=current, count=current, start=5, end=860);
        	filter.seqs(fasta=current, vertical=T, trump=.);
        	unique.seqs(fasta=current, count=current);
        	pre.cluster(fasta=current, count=current, diffs=2);
        	chimera.uchime(fasta=current, count=current, dereplicate=T);
        	remove.seqs(fasta=current, sranos=current);
        	classify.seqs(fasta=current, count=current, reference=data/references/trainset14_032015.pds.fasta, taxonomy=data/references/trainset14_032015.pds.tax, cutoff=80);
        	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);" '

rule calc_seq_dists:
    input:
        'data/raw/{sample}/{sample}.fasta'
    output:
        'data/processed/{sample}/{sample}.dist'
    params:
        mothur=mothur_bin,
        output_dir='data/processed/{sample}/'
    benchmark:
        'benchmarks/{sample}/calc_seq_dists.log'
    log:
        "logfiles/{sample}/calc_seq_dists.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.dir(output={params.output_dir}); dist.seqs(fasta={input[0]}, cutoff=0.03)"'
