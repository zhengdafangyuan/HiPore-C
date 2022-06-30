#!/bin/bash
SampleDir=$1
reffa=$2
fq=$3
threads=$4

sID=`basename $fq .fastq.gz`
Rawdir="${SampleDir}/Rawdata"
Mapdir="${SampleDir}/Mapping"
tempdir="${SampleDir}/tmpdir/${sID}"
mkdir -p $Mapdir
mkdir -p $tempdir

cd $tempdir

seqkit fx2tab ${Rawdir}/${fq} -i -n -l > "${sID}.fai" 
awk '{print $1}' "${sID}.fai" | sort > "${tempdir}/readslist.txt"

ngmlr -t $threads --subread-length 256 --max-segments 5 -r $reffa -q ${Rawdir}/${fq} -x ont -o "${tempdir}/reads_nglmr.sam" 
paftools.js sam2paf "${tempdir}/reads_nglmr.sam" > "${tempdir}/reads_nglmr.paf"  

# Unmapped Reads were realigned using minimap2
mappaf="${tempdir}/reads_nglmr.paf"
awk '{print $1}' ${mappaf} | sort | uniq > "${tempdir}/map_readslist.txt"
comm "${tempdir}/readslist.txt" "${tempdir}/map_readslist.txt" -3 | awk '{print $1}' > "${tempdir}/supp_readslist.txt"
seqtk subseq ${Rawdir}/${fq} "${tempdir}/supp_readslist.txt" > "${tempdir}/supp.fq"

suppfq="${tempdir}/supp.fq"
supppaf="${tempdir}/reads_supp.paf"
minimap2 -x map-ont -B 3 -O 2 -E 5 -k13 -t ${threads} -c $reffa "${suppfq}" > ${supppaf}

# Merge paf
awk -v OFS="\t" '($12>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$17}' $mappaf  > ${Mapdir}/${sID}.reads_map.paf
awk -v OFS="\t" '($12>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$22}' $supppaf >> ${Mapdir}/${sID}.reads_map.paf
