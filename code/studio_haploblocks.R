# Illumina WGS phased VCF
wgs <- read.delim('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/phased_wgs.gtf',as.is = T,header = F,stringsAsFactors = F)

# Nanopore phased VCF
nano <- read.delim('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/phased_Hg19_Nanopore.sort.gtf',as.is = T,header = F,stringsAsFactors = F)

# length haplotype blocks
summary(nano$V5-nano$V4)
summary(wgs$V5-wgs$V4)

if(F){
  
  # BED gene model
  bed <- read.delim('/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.bed',as.is = T,header = T,stringsAsFactors = F,check.names = F)
  summary(bed$end-bed$start)
  
  # filter BED gene model, keep only genes having a spanning haploblocks
  genes_wo_blocks <- read.delim('outs/tmp_analysis/genes_without_blocks.bed',as.is = T,header = F,stringsAsFactors = F)
  
  length(unique(genes_wo_blocks$V5))/length(unique(bed$GeneID)) # fraction of genes without haplotype blocks
  
  bed <- bed[which(!bed$GeneID %in% unique(genes_wo_blocks$V5)),] # filtered genes
  write.table(bed,file = "outs/tmp_analysis/humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.FILT.bed",col.names = T,sep = '\t',row.names = F,quote = F)
  
}

# intersect data
#intersectBed -a humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.FILT.bed -b phased_Hg19_Nanopore.sort.noChr.bed -wao > genes_intersect_haploblocks.bed
#intersectBed -a phased_Hg19_Nanopore.sort.noChr.bed -b phased_Hg19_Nanopore.sort.noChr.vcf -wa -wb > haploblocks_intersect_snps.bed

# process intersect data

genes_blocks <- read.delim('outs/tmp_analysis/genes_intersect_haploblocks.bed',as.is = T,header = F,stringsAsFactors = F,check.names = F)
genes_blocks <- genes_blocks[which( genes_blocks$V10 > 0 ),]
  
summary(genes_blocks$V10,100)

# check genes spanned by multiple haploblocks  
tocheck <- unique(genes_blocks$V5[which(duplicated(genes_blocks$V5))])

length(tocheck)/length(unique(genes_blocks$V5)) # fraction of genes spanned by multiple haploblocks

TabGenesBlocks <- genes_blocks[which(!genes_blocks$V5 %in% tocheck),]
TabGenesBlocks$haploblock <- paste(TabGenesBlocks$V6,TabGenesBlocks$V7,TabGenesBlocks$V8,sep = ':')

# phased SNPs per intersected gene and haplotypeblock 
gbs <- read.delim('outs/tmp_analysis/genes_intersect_haploblocks.bed',as.is = T,header = F,stringsAsFactors = F,check.names = F)

for(id in tocheck){
  message(id)
  count <- gbs[which(gbs$V5 == id),]
  count$haploblock <- paste(count$V6,count$V7,count$V8,sep = ':')
  
  this <- genes_blocks[which(genes_blocks$V5 == id),]
  this$haploblock <- paste(this$V6,this$V7,this$V8,sep = ':')
  
  block_to_add <- this[which(this$haploblock == count$haploblock[which.max(count$V10)]),]
  
  TabGenesBlocks <- rbind(TabGenesBlocks, block_to_add)
}


which(duplicated(TabGenesBlocks$V5)) # check for duplicates 

TabGenesBlocks <- TabGenesBlocks[order(TabGenesBlocks[,1],TabGenesBlocks[,2],TabGenesBlocks[,3]),]

# link haploblocks to snps
