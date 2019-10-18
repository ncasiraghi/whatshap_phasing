single_cells_bed <- list.files(path = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/tmp.parallel",pattern = 'cellid_',full.names = T)
# single_cells_bed <- list.files(path = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_GTseq/dna/cnvkit/wgs_normal_as_ctrl_bin_1Mb/cns_bed",pattern = 'cellid_',full.names = T)

haplotype_blocks <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.bed"

phased_snps <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf"

setwd('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks_10xCNV_data/9_clusters_362cells/')
# setwd('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks_GTseq_data/')

# intersect haplotype blocks and phased snps
cmd <- paste('intersectBed -a',haplotype_blocks,'-b',phased_snps,'-wa -wb > haplotype_blocks_with_phasedSNPs.bed')
system(cmd)

cmd <- paste('intersectBed -a',haplotype_blocks,'-b',phased_snps,'-c > haplotype_blocks_with_phasedSNPs_count.bed')
system(cmd)

x <- read.delim('haplotype_blocks_with_phasedSNPs_count.bed',stringsAsFactors = F,header = F)
summary(x$V5)
summary(x$V3-x$V2)

# n. of phased SNPs per haplotype block

# cluster cells based on Hana's clustering 
# scdna_clusters <- read.csv("/icgc/dkfzlsdf/analysis/B260/users/v390v/cnv_inference/data/raw/first_sample/scDNA/evo_dist_9_clustering.csv")
# scdna_clusters <- scdna_clusters[,2:4]

scdna_clusters <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV/STP/STP-Nuclei/processed_data/cluster9_cell_ids.txt")

scdna_clusters <- scdna_clusters[,2:3]
colnames(scdna_clusters) <- c("CELL_ID","LABEL")
barplot(table(scdna_clusters$LABEL),xlab = 'clusters (based on scDNA)',ylab = 'n. of cells') # n. cells per cluster. cluster 8 is the normal-cells cluster

scdna_bed <- data.frame(BED=single_cells_bed,CELL_ID=NA,stringsAsFactors = F)
scdna_bed$CELL_ID <- as.numeric(gsub(basename(scdna_bed$BED),pattern = 'cellid_|\\.bed$',replacement = ''))

scdna_bed <- merge(scdna_bed,scdna_clusters,by = 'CELL_ID',all.x = TRUE)
scdna_bed <- scdna_bed[which(!is.na(scdna_bed$LABEL)),]

clusters <- unique(scdna_bed$LABEL)

for(k in clusters){
  message(k)
  multiple_beds <- paste(scdna_bed$BED[which(scdna_bed$LABEL == k)],collapse = ' ')
  
  # step 1 : intersect all CN profiles from cells belonging to the same cluster
  step1 <- paste0('step1_',k,'.bed')
  cmd <- paste('bedops --partition',multiple_beds,'>',step1)
  system(cmd)
  
  # step 2 : intersect step1 with all BEDs to assign CN to each genomic segment 
  step2 <- paste0('step2_',k,'.bed')
  cmd <- paste('intersectBed -a',step1,'-b',multiple_beds,'-wa -wb | cut -f1,2,3,8 >',step2)
  system(cmd)
  
  # step 3 : 
  step3 <- paste0('step3_',k,'.bed')
  cmd <- paste('mergeBed -i',step2,'-d -1 -c 4 -o count,mean,median,collapse >',step3)
  system(cmd)
  
  m <- read.delim(file = step3,header = FALSE,as.is = TRUE,stringsAsFactors = FALSE)
  
  cn <- c()
  for(idx in seq_len(nrow(m))){
    cn <- c(cn, as.numeric(names(which.max(table(as.numeric(unlist(strsplit(m[idx,7],split = ','))))))))
  }
  
  m <- cbind(m,cn)
  
  # update step3
  write.table(m,file = step3,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  # step 4 : intersect haplotype blocks with step3 
  step4 <- paste0('step4_',k,'.bed')
  cmd <- paste('intersectBed -a',haplotype_blocks,'-b',step3,'-wa -wb >',step4)
  system(cmd)
  
  m <- read.delim(file = step4,header = FALSE,as.is = TRUE,stringsAsFactors = FALSE)
  
  m$haploblock <- paste(m[,1],m[,2],m[,3],sep = ":")
  
  hb <- unique(m$haploblock)
  
  hb_to_exclude <- c()
  tab <- c()
  for(x in hb){
    check <- length(unique(m[which(m$haploblock == x),12]))
    if( check > 1 ){
      hb_to_exclude <- c(hb_to_exclude, x)
    } else {
      tab <- rbind(tab, unique(m[which(m$haploblock == x),c(1:3,12)]))
    }
  }
  
  length(hb_to_exclude)/length(hb) # fraction of hb with multiple spanning CNs 
  
  m <- m[which(!m$haploblock %in% hb_to_exclude),c()]
  
  # step 5 : write haplo type blocks with associated copy number
  step5 <- paste0('step5_',k,'.bed')
  write.table(tab,file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  # step 6 : add phased SNPs to each block
  step6 <- paste0('step6_',k,'.bed')
  cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step6)
  system(cmd)
  
}
