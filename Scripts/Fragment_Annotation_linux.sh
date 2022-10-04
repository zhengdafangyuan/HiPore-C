#!/bin/bash
annopy=$1
SampleDir=$2
reffa=$3
fq=$4
gvdf=$5
ReAlignBash=$6
ChromCheckBash=$7

sID=`basename $fq .fastq.gz`
Rawdir="${SampleDir}/Rawdata"
Mapdir="${SampleDir}/Mapping"
tempdir="${SampleDir}/tmpdir/${sID}"
Exportdir="${SampleDir}/vdFAnnotation"
mkdir -p $Mapdir
mkdir -p $tempdir
mkdir -p ${Exportdir}

cd $tempdir
ln -s ${Rawdir}/${fq} ${tempdir}/${fq}
seqkit fx2tab ${tempdir}/${fq} -i -n -l > "${sID}.fai"

python $annopy "${tempdir}/${fq}" "${Mapdir}/${sID}.reads_map.paf" "${tempdir}" "${reffa}" "${gvdf}" "${ReAlignBash}" "${ChromCheckBash}" "${sID}.fai"
mv ${tempdir}/Align_Fragment_RvdF.csv "${Exportdir}/${sID}.Align_Fragment_RvdF.csv"
echo "$sID annotation finished!"
rm -rf ${tempdir}
