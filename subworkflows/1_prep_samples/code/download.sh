source /etc/profile.d/http_proxy.sh  # required for internet on the Great Lakes cluster
sra=$1
outdir=$2

prefetch $sra
fasterq-dump --split-files $sra -O $outdir
gzip ${outdir}/*.fastq