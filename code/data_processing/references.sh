#!/bin/bash

#Adapted from: https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017/blob/master/Makefile

REFS=data/references

wget -N -P $REFS/ http://www.mothur.org/w/images/b/be/Silva.nr_v123.tgz
tar xvzf $REFS/Silva.nr_v123.tgz -C $REFS/;
mothur "#get.lineage(fasta=$REFS/silva.nr_v123.align, taxonomy=$REFS/silva.nr_v123.tax, taxon=Bacteria-Archaea)"
mv $REFS/silva.nr_v123.pick.align $REFS/silva.bact_archaea.align
mv $REFS/silva.nr_v123.pick.tax $REFS/silva.bact_archaea.tax
rm $REFS/README.Rmd $REFS/README.html
rm $REFS/?ilva.nr_v123.*

mothur "#get.lineage(fasta=$REFS/silva.bact_archaea.align, taxonomy=$REFS/silva.bact_archaea.tax, taxon=Bacteria)"
mv $REFS/silva.bact_archaea.pick.align $REFS/silva.bacteria.align
mv $REFS/silva.bact_archaea.pick.tax $REFS/silva.bacteria.tax


mothur "#pcr.seqs(fasta=$REFS/silva.bacteria.align, start=13862, end=23445, keepdots=F, processors=8);degap.seqs();unique.seqs()"
cut -f 1 $REFS/silva.bacteria.pcr.ng.names > $REFS/silva.bacteria.pcr.ng.accnos
mothur "#get.seqs(fasta=$REFS/silva.bacteria.pcr.align, accnos=$REFS/silva.bacteria.pcr.ng.accnos);screen.seqs(minlength=240, maxlength=275, maxambig=0, maxhomop=8, processors=8); filter.seqs(vertical=T)"
mv $REFS/silva.bacteria.pcr.pick.good.filter.fasta $REFS/silva.v4.align
grep "^>" $REFS/silva.v4.align | cut -c 2- > $REFS/silva.v4.accnos
mothur "#get.seqs(taxonomy=$REFS/silva.bacteria.tax, accnos=$REFS/silva.v4.accnos)"
mv data/references/silva.bacteria.pick.tax data/references/silva.v4.tax
rm $REFS/silva.bacteria.pcr.*
rm $REFS/silva.filter


mkdir -p $REFS/rdp
wget -N -P $REFS/ http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
tar xvzf $REFS/Trainset14_032015.pds.tgz -C $REFS/rdp
mv $REFS/rdp/trainset14_032015.pds/trainset14_032015.* $REFS
rm -rf $REFS/rdp $REFS/Trainset*
