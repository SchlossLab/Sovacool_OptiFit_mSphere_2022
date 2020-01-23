source /etc/profile.d/http_proxy.sh  # required for internet on the Great Lakes cluster
list_file=$1
outdir=$2
prefetch --option-file $list_file
for SRA in $(cat $list_file)
do
    fasterq-dump --split-files $SRA -O $outdir
done
gzip ${outdir}/*