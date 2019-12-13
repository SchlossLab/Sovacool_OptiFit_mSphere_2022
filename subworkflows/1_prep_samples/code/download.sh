source /etc/profile.d/http_proxy.sh  # required for internet on the Great Lakes cluster
list_file=$1
outdir=$2
if [[ $outdir == *"mouse"* ]]; then
    wget -N -P $outdir http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar
    tar -xvf $outdir/StabilityNoMetaG.tar -C $outdir/
    rm $outdir/StabilityNoMetaG.tar
else
    for SRA in $(cat $list_file)
    do
        fasterq-dump --split-files $SRA -O $outdir
    done
    gzip ${outdir}/*
fi
