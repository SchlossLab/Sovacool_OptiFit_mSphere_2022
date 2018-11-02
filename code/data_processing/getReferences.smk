" Download the reference database and process with mothur "
import subprocess

configfile: 'code/data_processing/config_test.yaml'
output_dir = config['output_dir']
db_version = config['db_version']

rule all:
	input:
		expand("{output_dir}/silva.bacteria.{ext}", output_dir=output_dir, ext={'tax','align'}),
		expand("{output_dir}/trainset14_032015.pds.{ext}", output_dir=output_dir, ext={'tax','fasta'}),
		expand("{output_dir}/silva.v4.{ext}", output_dir=output_dir, ext={'accnos','align','tax'})

rule download_silva_db:
	output:
		"{{output_dir}}/Silva.nr_{version}.tgz".format(version=db_version)
	shell:
		'wget -N -P {{output_dir}} http://www.mothur.org/w/images/3/32/Silva.nr_{version}.tgz'.format(version=db_version)

rule unpack_silva_db:
	input:
		"{output_dir}/Silva.nr_{version}.tgz"
	output:
		"{output_dir}/silva.nr_{version}.align",
		"{output_dir}/silva.nr_{version}.tax"
	shell:
		"tar xvzf {{output_dir}}/Silva.nr_{version}.tgz -C {{output_dir}}/".format(version=db_version)

rule get_prok_lineage:
	input:
		fasta="{{output_dir}}/silva.nr_{version}.align".format(version=db_version),
		tax="{{output_dir}}/silva.nr_{version}.tax".format(version=db_version)
	output:
		"{output_dir}/silva.bact_archaea.align",
		"{output_dir}/silva.bact_archaea.tax"
	shell:
		"mothur '#get.lineage(fasta={{input.fasta}}, taxonomy={{input.tax}}, taxon=Bacteria-Archaea)' ; "
		"mv {{output_dir}}/silva.nr_{version}.pick.align {{output_dir}}/silva.bact_archaea.align ; "
		"mv {{output_dir}}/silva.nr_{version}.pick.tax {{output_dir}}/silva.bact_archaea.tax".format(version=db_version)

rule get_bact_lineage:
	input:
		fasta='{output_dir}/silva.bact_archaea.align',
		tax='{output_dir}/silva.bact_archaea.tax'
	output:
		'{output_dir}/silva.bact_archaea.pick.align',
		'{output_dir}/silva.bact_archaea.pick.tax'
	shell:
		"mothur '#get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria)'"

rule rename_bact:
	input:
		"{output_dir}/silva.bact_archaea.pick.{ext}"
	output:
		"{output_dir}/silva.bacteria.{ext}"
	shell:
		"mv {input} {output}"

rule pcr_seqs:
	input:
		"{output_dir}/silva.bacteria.align"
	output:
		"{output_dir}/silva.bacteria.pcr.ng.names",
		"{output_dir}/silva.bacteria.pcr.align"
	shell:
		'mothur "#pcr.seqs(fasta={input}, start=13862, end=23445, keepdots=F);degap.seqs();unique.seqs()"'

rule get_pcr_accession_numbers:
	input:
		"{output_dir}/silva.bacteria.pcr.ng.names"
	output:
		"{output_dir}/silva.bacteria.pcr.ng.accnos"
	shell:
		'cut -f 1 {input} > {output}'

rule get_fasta_seqs:
	input:
		fasta="{output_dir}/silva.bacteria.pcr.align",
		accnos="{output_dir}/silva.bacteria.pcr.ng.accnos"
	output:
		"{output_dir}/silva.bacteria.pcr.pick.good.filter.fasta"
	shell:
		'mothur "#get.seqs(fasta={input.fasta}, accnos={input.accnos});screen.seqs(minlength=240, maxlength=275, maxambig=0, maxhomop=8); filter.seqs(vertical=T)"'

rule rename_pcr_fasta:
	input:
		"{output_dir}/silva.bacteria.pcr.pick.good.filter.fasta"
	output:
		"{output_dir}/silva.v4.align"
	shell:
		'mv {input} {output}'

rule get_filtered_accession_numbers:
	input:
		"{output_dir}/silva.v4.align"
	output:
		"{output_dir}/silva.v4.accnos"
	shell:
		'grep "^>" {input} | cut -c 2- > {output}'

rule get_taxon_seqs:
	input:
		tax="{output_dir}/silva.bacteria.tax",
		accnos="{output_dir}/silva.v4.accnos"
	output:
		"{output_dir}/silva.bacteria.pick.tax"
	shell:
		'mothur "#get.seqs(taxonomy={input.tax}, accnos={input.accnos})"'

rule rename_taxon:
	input:
		"{output_dir}/silva.bacteria.pick.tax"
	output:
		"{output_dir}/silva.v4.tax"
	shell:
		'mv {input} {output}'

rule download_ribosomal_db:
	output:
		"{output_dir}/rdp/trainset14_032015.pds/trainset14_032015.pds.tax",
		"{output_dir}/rdp/trainset14_032015.pds/trainset14_032015.pds.fasta"
	shell:
		"wget -N -P {output_dir} http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz ; "
		"tar xvzf {output_dir}/Trainset14_032015.pds.tgz -C {output_dir}/rdp ; "

rule organize_ribosomal_db:
	input:
		expand("{{output_dir}}/rdp/trainset14_032015.pds/trainset14_032015.pds.{ext}", ext={'tax', 'fasta'})
	output:
		expand("{{output_dir}}/trainset14_032015.pds.{ext}", ext={'tax', 'fasta'})
	wildcard_constraints:
		output_dir="\w+"
	shell:
		"mv {output_dir}/rdp/trainset14_032015.pds/trainset14_032015.pds.* {output_dir}/ "
		#"rm -rf $REFS/rdp $REFS/Trainset*"
