# module load bedtools/2.24.0

library(data.table)
library(parallel)

# INPUTS scDNA 10x
if(FALSE){
  single_cells_bed <- list.files(path = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_10x/dna/sample1.genemodel.transcript/tmp.parallel",pattern = 'cellid_',full.names = T)
  cell_ids <- readLines(con = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks_10xCNV_data/scDNA_362cells_IDs.txt")
  outdir <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks_10xCNV_data/single_cells_362cells"
}

# INPUTS scDNA G&T
if(TRUE){
  single_cells_bed <- list.files(path = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_cn_GTseq/dna/cnvkit/wgs_normal_as_ctrl_bin_20kb_as_CellRanger/cns_bed/",pattern = '\\.nochr.bed$',full.names = T)
  cell_ids <- readLines(con = "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks_GTseq_data/scDNA_GTseq_89cells_IDs.txt")
  outdir <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/copy_number_blocks_GTseq_data/single_cells_89cells"
}

haplotype_blocks <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.3cols.bed"

phased_snps <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/tmp_analysis/phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf"

mc.cores <- 20

## RUN
setwd(outdir)

single_cells_bed <- single_cells_bed[which(basename(single_cells_bed) %in% paste0(cell_ids,".bed"))]

GetAnnotedSNPs <- function(i,single_cells_bed,haplotype_blocks,phased_snps){
  bed_file <- single_cells_bed[i]
  message(basename(bed_file))
  CELL_ID <- gsub(basename(bed_file),pattern = "\\.bed$",replacement = "")
  
  # step 4 : intersect haplotype blocks with step3 
  step4 <- paste0('step4_',CELL_ID,'.bed')
  cmd <- paste('intersectBed -a',haplotype_blocks,'-b',bed_file,'-wa -wb >',step4)
  system(cmd)
  
  m <- read.delim(file = step4,header = FALSE,as.is = TRUE,stringsAsFactors = FALSE)
  m <- m[,1:7]
  m$haploblock <- paste(m[,1],m[,2],m[,3],sep = ":")
  
  hb <- unique(m$haploblock)
  
  hb_to_exclude <- c()
  tab <- c()
  for(x in hb){
    check <- length(unique(m[which(m$haploblock == x),7]))
    if( check > 1 ){
      hb_to_exclude <- c(hb_to_exclude, x)
    } else {
      tab <- rbind(tab, unique(m[which(m$haploblock == x),c(1:3,7)]))
    }
  }
  
  length(hb_to_exclude)/length(hb) # fraction of hb with multiple spanning CNs 
  
  # step 5 : write haplo type blocks with associated copy number
  step5 <- paste0('step5_',CELL_ID,'.bed')
  write.table(tab,file = step5,col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  # step 7 : add phased SNPs to each block
  step7 <- paste0('step7_',CELL_ID,'.bed')
  cmd <- paste('intersectBed -a',step5,'-b',phased_snps,'-wa -wb >',step7)
  system(cmd)
  file.remove(step4,step5)
  
  # step 8 : reduce table columns and wrinte final output
  step8 <- paste0('step8_',CELL_ID,'.bed')
  cmd <- paste('cut -f 1-6,8,9,14',step7,'>',step8)
  system(cmd)
  file.remove(step7)
  
  m <- fread(file = step8,data.table = FALSE)
  header <- c("HB_CHROM","HB_START","HB_END","COPY_NUMBER","SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO")
  colnames(m) <- header
  
  reduce_info <- function(info){
    return(unlist(strsplit(info,split = ":"))[1])
  }
  
  m$SNP_PHASE_INFO <- unlist(lapply(X = m$SNP_PHASE_INFO,FUN = reduce_info))
  
  m <- m[,c("SNP_CHROM","SNP_POS","SNP_REF","SNP_ALT","SNP_PHASE_INFO","HB_CHROM","HB_START","HB_END","COPY_NUMBER")]
  write.table(m,file = step8,col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  
}

mclapply(seq(single_cells_bed),GetAnnotedSNPs,single_cells_bed=single_cells_bed,haplotype_blocks=haplotype_blocks,phased_snps=phased_snps,mc.cores = mc.cores)
