#!/bin/bash
rawfq=$1
chrname=$2
subreadbed=$3
paffile=$4
workdir=$5
rawfq=$1

# seperated chromosome fastas
chrrefDir="~/GCA_000001405.15_GRCh38_no_alt_analysis_set"
chrfa="${chrrefDir}/${chrname}.fa"

mkdir -p $workdir

cd ${workdir}
seqtk subseq $rawfq $subreadbed >  ${workdir}/${chrname}.subreads.fq

minimap2 -x map-ont -k13 -t 5 $chrfa ${workdir}/${chrname}.subreads.fq > ${paffile}
rm ${workdir}/${chrname}.subreads.fq
