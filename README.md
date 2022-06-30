# HiPore-C
We developed a protocol of in situ high throughput multi-way contact long read Pore-C sequencing (in situ HiPore-C), a strategy that integrated multi-fragment ligates preparation with third-generation sequencing technology. Compared to the reported Pore-C method, HiPore-C can yield more chromatin interactions than traditional Hi-C and Pore-C at the same cost using simple procedures.And base on the high-order and allele-specific feature of HiPore-C, we have globally characterized single-allele topologies with unprecedented depth to reveal elusive genome folding principles.

This is the code used to make a pipeline for analysing in situ high throughput multi-way contact long read Pore-C sequencing (HiPore-C) data. In this work, We use the human hg38 genome as reference, and we obtained the public Hi-C, Chip-seq, DNase-seq and RNA-seq datasets of GM12878 and K562 cell lines from 4DN or ENDCODE database. 

# Software
This pipeline has a number of dependencies including the following:

Guppy (4.5.3);  
Megalodon (2.3.4);  
Nanoplot (1.30.1);  
ngmlr (0.2.7);  
minimap2 (2.17-r941);  
samtools (1.10);  
gffcompare (version 0.11.6);  
Pore-C-Snakemake pipeline (0.3.0);  
FIMO (5.3.3);   
cooltools (0.5.1);  
cooler (0.8.6);  
HiCrep (1.2.0);  
Juicer(1.22.01);  
Juicebox(v2.10.01);  
HiCExplorer (3.6);  
HiGlass  
FAN-C 


# Basecalling and Methylation calling

In steps, Guppy and Megalodon software were used. And due to the large amount of Nanopore data, it is recommended to generate a corresponding fastq or modification callling bam result for each multfast5 file (typically containing 4000 reads) in order to facilitate the subsequent analysis. The input is fast5 file and the output is fastq and mod bam file.

``` 
Conf_file="~/ont-guppy/data/dna_r9.4.1_450bps_hac_prom.cfg"

guppy_basecaller \
-i ${Fast5_file_dir} \
-s ${Exportdir} \
-c ${Conf_file}  \
-x cuda:0 \
-r -q 0 \
--min_qscore 7 \
--compress_fastq
``` 

``` 
megalodon \
    ${Fast5_file_dir} \
    --guppy-server-path "~/ont-guppy/bin/guppy_basecall_server" \
    --guppy-params "-d ~/rerio/basecall_models/" \
    --guppy-config res_dna_r941_prom_modbases_5mC_v001.cfg \
    --outputs mod_basecalls \
    --output-directory "${Exportdir}" \
    --mod-motif m CG 0 \
    --devices all --processes 48 --overwrite
``` 

# Mapping and Fragment Annotation
In this step, to obtain more accurate alignment results for high-order reads, and to improve the effective usage of HiPore-C data, we introduced the NGMLR and Minimap2 for three-generation sequencing to construct the HiPore-C alignment pipeline. NGMLR algorithm has advantages for sequence breakpoint identification and Minimap2 algorithm has advantages for long noisy sequence reads alignment. Consistent with the above, we recommend to perform alignment and annotation for each single fastq file (generated from a 4000 reads multifast5 file).

Use the ./Scripts/Alignment.sh  and ./Scripts/Fragment_Annotation.sh script to perform this analysis. The fastq data needs to be placed in the ```${Sampledir}/Rawdata``` directory, the generated alignment result files will be placed to ```${Sampledir}/Mapping``` , and the annotation files will be placed to ```${Sampledir}/vdFAnnotation``` directory.  A python script for annotating fragments ```annopy="./Scripts/Read_Fragment_Annotation.py"```, a table of in-silicon genomic restriction enzyme digested fragments ```genome_digest_frag="./Scripts/DpnII_GRCh38.vd.fragments.csv"``` , and a shell script for re-alignment ```ReAlignBash="./Scripts/minimap2_subreads_remapping.sh"``` and alignemtn validation ```ChromCheckBash="./Scripts/multichrom_check.sh"``` are also used in this step of the analysis.

To run the pipeline type:

``` 
base ./Scripts/Alignment.sh $Sampledir $refgenome $fastq $threads
base ./Scripts/Fragment_Annotation.sh $annopy $Sampledir $refgenome $fastq $genome_digest_frag  $ReAlignBash $ChromCheckBash

``` 


# Generate pairwise contact matrix
To compare with the Hi-C data and the subsequent visualization, the multiway contacts of HiPore-C reads were decomposed into pairwise contacts using ```Generate_Contact_juiceMatrix.py```， and then  mcool and hic format file were generated by cooler and Juicer tools, repectively.

To generate pairwise contact matrix, merge all ```Align_Fragment_RvdF.csv``` files from above steps, and type:
``` 
inputfile="Merge_Align_Fragment_RvdF.csv"
juice_matrix="juice_matrix.txt"
python ./Scripts/Generate_Contact_juiceMatrix.py ${inputpaf} ${contact_matrix}
``` 


To generate .cool, .mcool and .hic files type:
``` 
DataDir="${SampleDir}/juiceMatrix"
PAIRS_FILE="${DataDir}/juice_matrix.txt"
PAIRS_FILE_Sort="${DataDir}/juice_matrix_sort.txt"

CHROMSIZES_FILE="~/hg38.chromosomes.size"
BINSIZE=1000
COOL_FILE="${DataDir}/output.cool"
HIC_FILE="${DataDir}/output.hic"

# sorted
awk '{if ($3 > $7) { print $1, $6, $7, $8, $9, $2, $3, $4, $5, $11, $10} else {print $0} }' ${PAIRS_FILE} | sort -k3,3d -k7,7d | tr " " "\t"  >  ${PAIRS_FILE_Sort}

# contact matrix to cool and mcool
cooler cload pairs -c1 3 -p1 4 -c2 7 -p2 8 $CHROMSIZES_FILE:$BINSIZE $PAIRS_FILE_Sort $COOL_FILE && \
cooler zoomify -p 20 -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000 --balance ${COOL_FILE} &

# juice_matrix to hic
java -Xmx20g -jar ~/juicer_tools_1.22.01.jar pre --threads 12 ${PAIRS_FILE_Sort} ${HIC_FILE} hg38 &
``` 


# High order chromatin structure analysis

For the validation of chromation structure, we analyszed the compartment eigenvector scores, insulation scores and performed APA analysis.
``` 
DataDir="${SampleDir}/juiceMatrix"
mCOOL_FILE="${DataDir}/output.mcool"
mCOOL_FILE_10kb="${mCOOL_FILE}::resolutions/10000"
mCOOL_FILE_100kb="${mCOOL_FILE}::resolutions/100000"

# Eigenvector and Insulation scores
cooltools eigs-cis --bigwig -o Cis_Compartment.100k $mCOOL_FILE_100kb &
cooltools insulation -o Insulation.tsv --window-pixels --append-raw-scores --bigwig $mCOOL_FILE_10kb 1 2 5 10 &

# APA
loopfile="~/Rao2014_Loops/GSE63525_GM12878_hg38_looplist.txt"
hicfile="${SampleDir}/juiceMatrix/out.hic"
Exportdir="${SampleDir}/juiceMatrix/APA"
mkdir -p ${Exportdir}
java -Xmx20g -jar ~/juicer_tools_1.22.01.jar apa --threads 20 -r 10000 -w 10 -u ${hicfile} ${loopfile} ${Exportdir}
``` 
