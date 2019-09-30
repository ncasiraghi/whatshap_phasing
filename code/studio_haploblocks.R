# Illumina WGS phased VCF
wgs <- read.delim('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/phased_wgs.gtf',as.is = T,header = F,stringsAsFactors = F)

# Nanopore phased VCF
nano <- read.delim('/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/whatshap_phasing/outs/phased_Hg19_Nanopore.sort.gtf',as.is = T,header = F,stringsAsFactors = F)

# length haplotype blocks
summary(nano$V5-nano$V4)
summary(wgs$V5-wgs$V4)


# process data

# intersectBed -a phased_Hg19_Nanopore.sort.noChr.bed -b /icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13_TranscriptStartEnd.sort.merged.unique.bed -wb > blocks_int_genes.bed
blocks_int_genes <- read.delim('outs/tmp_analysis/blocks_int_genes.bed',as.is = T,header = F,stringsAsFactors = F,check.names = F)
summary(blocks_int_genes$V3-blocks_int_genes$V2)
blocks_int_genes <- blocks_int_genes[which(blocks_int_genes$V3-blocks_int_genes$V2 > 0),]
blocks_int_genes <- blocks_int_genes[order(blocks_int_genes[,1],blocks_int_genes[,2],blocks_int_genes[,3]),]

# count n. snps per haploblock
# intersectBed -a blocks_int_genes.bed -b phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf -c > blocks_int_genes_SNPs_count.bed
snps_count <- read.delim('outs/tmp_analysis/blocks_int_genes_SNPs_count.bed',as.is = T,header = F,stringsAsFactors = F,check.names = F)

# check genes spanned by multiple haploblocks  
tocheck <- unique(blocks_int_genes$V8[which(duplicated(blocks_int_genes$V8))])

length(tocheck)/length(unique(blocks_int_genes$V8)) # fraction of genes spanned by multiple haploblocks

TabBlocksGenes <- blocks_int_genes[which(!blocks_int_genes$V8 %in% tocheck),]

for(id in tocheck){
  message(id)
  m <- snps_count[which(snps_count$V8 == id),]
  m <- m[which.max(m$V10),]
  group <- paste(m$V1,m$V2,m$V3,sep = ":")
  all <- paste(blocks_int_genes$V1,blocks_int_genes$V2,blocks_int_genes$V3,sep = ":")
  TabBlocksGenes <- rbind(TabBlocksGenes,blocks_int_genes[which(all==group),])
}


which(duplicated(TabBlocksGenes$V8)) # sanity check for duplicates 

TabBlocksGenes <- TabBlocksGenes[order(TabBlocksGenes[,1],TabBlocksGenes[,2],TabBlocksGenes[,3]),]

write.table(x = TabBlocksGenes,file = "outs/tmp_analysis/TabBlocksGenes.bed",quote = F,col.names = F,row.names = F,sep = "\t")

# intersect TabBlocksGenes and phased SNPs
# intersectBed -a TabBlocksGenes.bed -b phased_Hg19_Nanopore.sort.noChr.OnlyPhased.vcf -wa -wb > TabBlocksGenes_with_phasedSNPs.bed

