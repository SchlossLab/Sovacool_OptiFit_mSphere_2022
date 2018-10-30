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
		"README.md",
		"silva.nr_{version}.align",
		"silva.nr_{version}.tax"
	shell:
		"tar xvzf {dir}/Silva.nr_{{version}}.tgz -C {dir}/".format(dir=output_dir)

'''
rule get_prokaryotic_lineage:
	input:
	output:
	shell:
		"mothur '#get.lineage(fasta={dir}/silva.nr_{version}.align, taxonomy={dir}/silva.nr_{version}.tax, taxon=Bacteria-Archaea)'".format(dir=output_dir, version=db_version)
'''
