" Download the reference database and process with mothur "
import subprocess

configfile: 'code/data_processing/config_test.yaml'
output_dir = config['output_dir']
db_version = config['db_version']

rule all:
	input:
		expand("{output_dir}/silva.bact_archeae.{ext}", output_dir=output_dir, ext={'tax', 'align'})

rule download:
	output:
		"{{output_dir}}/Silva.nr_{version}.tgz".format(version=db_version)
	shell:
		'wget -N -P {output_dir} http://www.mothur.org/w/images/3/32/Silva.nr_{version}.tgz'

rule untar:
	input:
		"{output_dir}/Silva.nr_{version}.tgz"
	output:
		"{output_dir}/silva.nr_{version}.align",
		"{output_dir}/silva.nr_{version}.tax"
	shell:
		"tar xvzf {output_dir}/Silva.nr_{version}.tgz -C {output_dir}/"

rule get_prok_lineage:
	input:
		align="{output_dir}/silva.nr_{version}.align",
		tax="{output_dir}/silva.nr_{version}.tax"
	output:
		"{output_dir}/silva.nr_{version}.pick.align",
		"{output_dir}/silva.nr_{version}.pick.tax"
	shell:
		"mothur '#get.lineage(fasta={input.align}, taxonomy={output_dir}/silva.nr_{version}.tax, taxon=Bacteria-Archaea)'"

rule rename:
	input:
		"{{output_dir}}/silva.nr_{version}.pick.{{ext}}".format(version=db_version)
	output:
		"{output_dir}/silva.bact_archeae.{ext}"
	shell:
		"mv {input} {output}"
