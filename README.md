# Assign copy number status and phased SNPs to Haplotype blocks.

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
The VCF output from WhatsHap.
