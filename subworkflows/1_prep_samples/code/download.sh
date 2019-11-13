list_file=$1
outdir=$2
if [[ $outdir == *"mice"* ]]; then
    wget -N -P $outdir http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar
    tar -xvf $outdir/StabilityNoMetaG.tar -C $outdir/
    rm $outdir/StabilityNoMetaG.tar
else
    for SRA in $(cat $list_file)
    do
        prefetch $SRA
        fasterq-dump --split-files $SRA -O $outdir
    done
    gzip ${outdir}/*
fi
