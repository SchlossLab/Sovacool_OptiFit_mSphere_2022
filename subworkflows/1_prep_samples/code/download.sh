list=$1
outdir=$2
if [[ $outdir == *"mice"* ]]; then
    wget -N -P $outdir http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar
    tar -xvf $outdir/StabilityNoMetaG.tar -C $outdir/
    rm $outdir/StabilityNoMetaG.tar
else
    for SRA in $(cat $list)
    do
        fastq-dump --split-files $SRA -O $outdir --gzip
    done
fi
