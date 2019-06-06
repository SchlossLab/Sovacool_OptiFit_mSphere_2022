output_dir = os.path.join(config['input_dir'], 'references')

rule rdp_targets:
    input:
        expand("{output_dir}/rdp/trainset14_032015.pds/trainset14_032015.pds.{ext}", output_dir=output_dir, ext={'tax', 'fasta'})

rule download_rdp_db:
    output:
        expand("data/references/rdp/trainset14_032015.pds/trainset14_032015.pds.{ext}",ext={'tax', 'fasta'})
    params:
        dir="data/references/rdp/",
        tar="Trainset14_032015.pds.tgz"
    shell:
        "wget -N -P {params.dir} http://www.mothur.org/w/images/8/88/{params.tar} ; "
        "tar xvzf {params.dir}{params.tar} -C {params.dir} ; "
