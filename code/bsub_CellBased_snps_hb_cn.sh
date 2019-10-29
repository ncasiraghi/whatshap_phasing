#!/bin/sh 
#BSUB -J scBED
#BSUB -q long
#BSUB -e scBED.log 
#BSUB -o scBED.txt 

module load R/3.6.0 bedtools/2.24.0

Rscript /icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/code/CellBased_snps_hb_cn.R