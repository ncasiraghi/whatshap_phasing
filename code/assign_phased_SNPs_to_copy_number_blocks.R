single_cells_bed <- list.files(path = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/tmp.parallel",pattern = 'cellid_',full.names = T)

haplotype_blocks <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.bed"

phased_snps <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf"

setwd('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks/')

# intersect sc copy number block with haplotype blocks and assign phased snps to each segment

if(FALSE){
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
}

# n. of phased SNPs per haplotype block

# cluster cells based on Hana's clustering 
scdna_clusters <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/cnv_inference/data/raw/first_sample/scDNA/evo_dist_9_clustering.csv")
scdna_clusters <- scdna_clusters[,2:4]

barplot(table(scdna_clusters$LABEL),xlab = 'clusters (based on scDNA)',ylab = 'n. of cells') # n. cells per cluster. cluster 8 is the normal-cells cluster

scdna_bed <- data.frame(BED=single_cells_bed,CELL_ID=NA,stringsAsFactors = F)
scdna_bed$CELL_ID <- as.numeric(gsub(basename(scdna_bed$BED),pattern = 'cellid_|\\.bed$',replacement = ''))

scdna_bed <- merge(scdna_bed,scdna_clusters,by = 'CELL_ID',all.x = TRUE)
scdna_bed <- scdna_bed[which(!is.na(scdna_bed$LABEL)),]

clusters <- unique(scdna_bed$LABEL)

for(k in clusters){
  message(k)
  multiple_beds <- paste(scdna_bed$BED[which(scdna_bed$LABEL == k)],collapse = ' ')
  bedops_out <- paste0('bedops_partition_cluster_',k,'.bed')
  cmd <- paste('bedops --partition',multiple_beds,'>',bedops_out)
  system(cmd)
  
  intersectbed_out <- paste0('intersectbed_out_',k,'.bed')
  cmd <- paste('intersectBed -a',bedops_out,'-b',multiple_beds,'-wa -wb >',intersectbed_out)
  system(cmd)
  
  intersectbed_out_reduced <- gsub(intersectbed_out,pattern = "\\.bed",replacement = '_reduced.bed')
  cmd <- paste('cut -f1,2,3,8',intersectbed_out,'>',intersectbed_out_reduced)
  system(cmd)
  
  mergebed_out <- paste0('mergebed_out_',k,'.bed')
  cmd <- paste('mergeBed -i',intersectbed_out_reduced,'-c 4 -o median >',mergebed_out)
  system(cmd)
  
  # intersect merged bed with haplotype-block 
  final <- paste0('final_',k,'.bed')
  cmd <- paste('intersectBed -a',haplotype_blocks,'-b',mergebed_out,'-wb >',final)
  system(cmd)
  
  final_reduced <- paste0('final_reduced_',k,'.bed')
  cmd <- paste('cut -f 1,2,3,8',final,'>',final_reduced)
  system(cmd)
  
  # add phased SNPs to each block
  tabfinal <- paste0('TabHaplotypeblock_with_phasedSNPs_',k,'.bed')
  cmd <- paste('intersectBed -a',final_reduced,'-b',phased_snps,'-wa -wb >',tabfinal)
  system(cmd)
  
}


