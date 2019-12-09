#!/bin/bash

#------------------------------------------
##SNPs qc
#wd=/gpfs23/data/coxvgi/zhoud2/data/gtex/geno/v6/
#plink --bfile ${wd}/maf0.05 --geno 0.05 --hwe 0.05 --make-bed --out ${wd}/v6p

##plink bed to dosage
#wd=/data/coxvgi/zhoud2/data/gtex/geno
#mkdir ${wd}/v6/dosage
#python2.7 /home/zhoud2/tmp/predixcan/convert_plink_to_dosage_backup.py -p plink -b ${wd}/v6/v6p -o ${wd}/v6/dosage/v6p

#generate genotype file for each gene
#ml GCC OpenMPI R
#Rscript /home/zhoud2/script/neanderthals/geno_eachgene.r ${SLURM_ARRAY_TASK_ID} #args=1-22

#generate samples with Neandthals haplotypes for each gene
#ml GCC OpenMPI R
#Rscript /home/zhoud2/script/neanderthals/sample_list.r

#gene expression prediction model training and evaluation
ml GCC OpenMPI R
#Rscript /home/zhoud2/script/neanderthals/en.r ${SLURM_ARRAY_TASK_ID}
Rscript /home/zhoud2/script/neanderthals/en_resample.r ${SLURM_ARRAY_TASK_ID}

#combine results
#sh /home/zhoud2/script/neanderthals/combine.sh

#plot results
#ml GCC OpenMPI R
#Rscript /home/zhoud2/script/neanderthals/plot.r


