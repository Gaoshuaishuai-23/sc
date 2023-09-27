#!/usr/bin/bash

#PBS -N hro5
#PBS -l walltime=24:00:00 
#PBS -q workq
#PBS -l mem=100GB
#PBS -l ncpus=32

memo=HRO-5
#DATA=/share2/pub/zhouyj/zhouyj/organoid/data/eye/analysis/${memo}
DADA2=/share2/pub/zhouyj/zhouyj/organoid/data/eye/analysis/${memo}/drived_data

#mkdir $DATA
#mkdir $DADA2

cd $DARA2

source /share/pub/zhouyj/anaconda3/bin/activate SC

#Rscript /share2/pub/zhouyj/zhouyj/organoid/data/eye/E_MTAB_12714/script/scATAC_5.r
#Rscript /share2/pub/zhouyj/zhouyj/organoid/data/eye/E_MTAB_12714/script/enrich_ATAC_5.r
#python /share2/pub/zhouyj/zhouyj/organoid/data/eye/E_MTAB_12714/script/ArchR_h5ad_5.py

