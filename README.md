## Assign copy number and phased SNPs to Haplotype blocks.

### Input 1. Copy number profiles
Each entry indicates the copy number status of a genomic segment.<br/>
The format is simply 4 columns, tab-delimited, without header.
```
column 1: chromosome
column 2: start
column 3: end
column 4: integer copy number 
```
### Input 2. Haplotype blocks
Each entry is a Haplotype Block.<br/>
The format is simply 3 columns, tab-delimited, without header.
```
column 1: chromosome
column 2: start
column 3: end
```
### Input 3. Phased SNPs
Each entry is a SNP with phased information such as 0|1 or 1|0.<br/>
The VCF output from WhatsHap.

### Input 4. Cluster assignment
Each entry is a cell with associated cluster label.<br/>
The format is simply 2 columns, tab-delimited, without header.
```
column 1: cell ID
column 2: cluster label
```
