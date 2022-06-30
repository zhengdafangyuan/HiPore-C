#!/bin/bash
workdir=$1
rawfq=$2
subreadbed=$3
paffile=$4
reffa=$5

mkdir -p $workdir
cd ${workdir}
seqtk subseq $rawfq $subreadbed >  ${workdir}/subreads.fq
minimap2 -x map-ont -B 3 -O 2 -E 5 -k13 -t 22 $reffa ${workdir}/subreads.fq > ${paffile}
rm ${workdir}/subreads.fq
