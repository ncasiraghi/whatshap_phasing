single_cells_bed <- list.files(path = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/tmp.parallel",pattern = 'cellid_',full.names = T)

haplotype_blocks <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.bed"

phased_snps <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf"

setwd('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks/')

# intersect sc copy number block with haplotype blocks and assign phased snps to each segment

for(sc in single_cells_bed){
  message(basename(sc))
  outbed <- gsub(basename(sc),pattern = '\\.bed$',replacement = '_intersect_haploblock.bed')
  cmd <- paste('intersectBed -a',sc,'-b',haplotype_blocks,'-wb >',outbed)
  system(cmd)
  
  haploblock_cn <- read.delim(outbed,as.is = T,header = F,stringsAsFactors = F,check.names = F)
  summary(haploblock_cn$V3-haploblock_cn$V2)
  
  haploblock_cn <- haploblock_cn[which(haploblock_cn$V3-haploblock_cn$V2 > 0),]
  haploblock_cn <- haploblock_cn[order(haploblock_cn[,1],haploblock_cn[,2],haploblock_cn[,3]),]
  
  write.table(haploblock_cn,file = outbed,col.names = F,row.names = F,quote = F,sep = "\t")
  
  finalbed <- gsub(basename(sc),pattern = '\\.bed$',replacement = '_haploblocks_with_phasedSNPs.bed')
  cmd <- paste('intersectBed -a',outbed,'-b',phased_snps,'-wa -wb >',finalbed)
  system(cmd)
  
}

# n. of phased SNPs per haplotype block

# cluster cells based on cn-phased-blocks
