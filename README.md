# HiPore-C
We developed a protocol of in situ high throughput multi-way contact long read Pore-C sequencing (in situ HiPore-C), a strategy that integrated multi-fragment ligates preparation with third-generation sequencing technology. With HiPore-C approach, we could explore higher-order chromatin interaction genome-widely.

# Introduction

in situ HiPore-C

# Software
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

# Basecalling
``` 
Conf_file="~/ont-guppy/data/dna_r9.4.1_450bps_hac_prom.cfg"

guppy_basecaller \
-i ${Fast5dir} \
-s ${Exportdir} \
-c ${Conf_file}  \
-x cuda:0 \
-r -q 0 \
--min_qscore 7 \
--compress_fastq
``` 

# Methylation calling
``` 
megalodon \
    ${Fast5dir} \
    --guppy-server-path "./ont-guppy/bin/guppy_basecall_server" \
    --guppy-params "-d ./rerio/basecall_models/" \
    --guppy-config res_dna_r941_prom_modbases_5mC_v001.cfg \
    --outputs mod_basecalls \
    --output-directory "${Exportdir}" \
    --mod-motif m CG 0 \
    --devices all --processes 48 --overwrite
``` 

# Mapping and Fragment Annotation

Aignment.sh

``` 
``` 

# Generate pairwise contact matrix

Generate_Contact_juiceMatrix.py
``` 
inputfile="Merge_Align_Fragment_RvdF.csv"
juice_matrix="juice_matrix.txt"
python ./Scripts/Generate_Contact_juiceMatrix_v1.1py ${inputpaf} ${contact_matrix}
``` 

# High order chromatin structure analysis
