" Download the reference database and process with mothur "
import subprocess

config: config.yaml
output_dir = config['output_dir']
db_version = config['db_version']

rule download:
	output:
		"{dir}/Silva.nr_{version}.tgz".format(dir=output_dir, version=db_version)
	shell:
		'wget -N -P {dir} http://www.mothur.org/w/images/3/32/Silva.nr_{version}.tgz'.format(dir=output_dir, version=db_version)

rule untar:
	input:
		"{dir}/Silva.nr_{{version}}.tgz".format(dir=output_dir)
	output:
		"{dir}/README.md",
		"{dir}/silva.nr_{version}.align",
		"{dir}/silva.nr_{version}.tax"
	shell:
		"tar xvzf {dir}/Silva.nr_{{version}}.tgz -C {dir}/".format(dir=output_dir)

rule get_prok_lineage:
	input:
		readme="{dir}/README.md"
		align="{dir}/silva.nr_{version}.align",
		tax="{dir}/silva.nr_{version}.tax"
	output:
		"{dir}/silva.nr_{version}.pick.align",
		"{dir}/silva.nr_{version}.pick.tax"
	shell:
		"mothur '#get.lineage(fasta={input.align}, taxonomy={dir}/silva.nr_{version}.tax, taxon=Bacteria-Archaea)'"

'''
rule rename:
	input:
		align="{dir}/silva.nr_{VERSION}.pick.align",
		tax="{dir}/silva.nr_{VERSION}.pick.tax"
	output:
		align="{dir}/silva.bact_archaea.align",
		tax="{dir}/silva.bact_archaea.tax"
	shell:
		"mv {input.align} {output.align} ; "
		"mv {input.tax} {output.tax}"
'''

rule rename:
	input:
		"{dir}/silva.nr_{version}.pick.{ext}"
	output:
		"{dir}/silva.bact_archeae.{ext}"
	shell:
		"mv {input} {output}"
