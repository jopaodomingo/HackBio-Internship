# Project 3 (RNA-Seq Plants)
## Introduction
In this project, analysis of the <i>Arabidopsis thaliana</i> leaf's response to UV stress in the vasculature was conducted using Bash and RStudio. The objectives include the following: 
1) To perform differential expression analysis in the vasculature of UV treated <i>Arabidopsis thaliana</i> tissues with Water (control) and UV-C light
2) To determine the top 100 differentially expressed genes of the two groups of <i>Arabidopsis thaliana</i> tissues
3) To determine the top 5 enriched pathways of the two groups of <i>Arabidopsis thaliana</i> tissues

## Methodology
### Dataset
The data set consisted of 6 samples that were obtained from SRA-Explorer. These samples consisted of 3 replicates that were exposed to water (SRR12808527, SRR12808528, SRR12808529) and UV-C light (SRR12808497, SRR12808498, SRR12808499). A series of pre-processing steps were conducted as follows 

### FastQC implementation
```
#!/bin/bash

#Quality Control
mkdir -p qc

for filename in clean_files/*.gz; do
        fastqc -o qc/ $filename
        gunzip "$filename"
done

#MultiQC Compilation
multiqc qc/ -o qc/
```
<img width="1506" height="813" alt="image" src="https://github.com/user-attachments/assets/9dd9de2b-0e41-40a0-8f12-6f84fc10f6d1" />
Issues in the per base sequence content, sequence duplication levels, and overrepresented sequences of the samples suggested the need for trimming via FastP.

### FastP implementation
```
#!/bin/bash 
#Trimming Reads
mkdir -p trim

for filename in clean_files/*.fastq; do
        base=$(basename "$filename" .fastq)
        fastp -i "$filename" -o "trim/${base}_trim.fastq" -h "trim/${base}_report.html"
done
```
After trimming, reference genome mapping was done using STAR. The reference genome was obtained from Ensembl. 
### Reference genome mapping using STAR
> nano star_map.sh
```
#!/bin/bash

#make directory for genome annotation
mkdir genome

#download reference genome
curl -L https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -o genome/a_thaliana.fa.gz
gunzip genome/a_thaliana.fa.gz

# create genome index directory and assemble genome
STAR --runMode genomeGenerate --genomeDir genome/genomeIndex --genomeFastaFiles genome/a_thaliana.fa

# create directories for output
mkdir -p mapped

# Genome mapping proper
for infile in trim/*.fastq ; do
        outfile=$(basename "$infile" .fastq)
        STAR --genomeDir genome/genomeIndex --readFilesIn $infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

# Create an IGV directory
mkdir IGV

# Index the bam files for IGV and create .bai file

for infile in trim/*.bam ; do
        samtools index -o IGV/ $infile
done

```
The tool featureCounts is then used to generate the data containing the upregulated and downregulated gene profile:
### featureCounts
> nano feature_count.sh
```
#!/bin/bash

#featurecounts and reads
mkdir counts

#download genome annotation
curl -L https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz -o genome/a_thaliana.gff3.gz
gunzip genome/a_thaliana.gff3.gz

#run featureCounts
featureCounts -O -t gene -g ID -a genome/a_thaliana.gff3 -o counts/counts.txt IGV/*.bam
```
The data is then uploaded to RStudio, where Differential Sequencing Analysis (DSeq2) is applied. 
```
#set working directory
setwd('/Users/Jopao/Desktop/Graduate School Applications/HackBio/Project 3')

#libraries
library(DESeq2) #For Differential Expression Analysis
library(pheatmap) #For generation of heat map

#count file; files must be found in wd
c_e_count <- read.delim('counts.txt', header = T)
c_e_meta <- read.delim('metadata.txt', header = T, stringsAsFactors = T)

#preview (recheck the first few lines)
head(c_e_count)
head(c_e_meta)

#keep to the important columns; only consider columns SRR12808497 to SRR12808529
raw_counts <- c_e_count[, c("SRR12808497", "SRR12808498", "SRR12808499", "SRR12808527", "SRR12808528", "SRR12808529")]
head(raw_counts)

#add rownames; change row numbers to the gene IDs
rownames(raw_counts) <- c_e_count$Geneid
head(raw_counts)

#create desqdataset; countData is the actual data; colData is the metadata for differential expression analysis 
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = c_e_meta,
                              design = ~stress_exposure)

#preview
dds
dds$sample
dds$stress_exposure

#Perform the differential expression analysis
dds <- DESeq(dds)

#Assign the final results to variable final_res
final_res <- results(dds)

#look at your result
head(final_res)

#we have a truncated data, let's see the distribution of p-values
plot(density(x = na.omit(final_res$pvalue)))

#ok let's look at our differentially expressed genes
plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$padj),#x is the magnitude of change, y is the distribution of p values 
     cex = 0.25,#size of characters
     pch = 19, #shape used 
     col = 'grey', #color
     ylim = c(0,20), #y-range 
     ylab = 'Adjusted P-Value', #label
     xlab = 'Log2 FC') #label
     
     #Generate absolute vertical lines and a horizontal line (h)
     abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)
     
     #where are the upregulated (increased genes in the male relative to female)
     upregulated <- subset(final_res, padj < 0.05 & log2FoldChange > 2) #log2FoldChange is change in magnitude
     points(upregulated$log2FoldChange,
            y = -log10(upregulated$padj), 
            cex = 0.35,
            pch = 19,
            col = 'salmon')
     
     #where are the downregulated
     downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -2)
     points(downregulated$log2FoldChange,
            y = -log10(downregulated$padj), 
            cex = 0.35,
            pch = 19,
            col = 'lightblue')
     
     mtext('Volcano plot of gene regulation in UV-C stressed Arabidopsis thaliana')
     #we can merge the two to do a clean and less memory efficient heatmap
     degs <- rbind(raw_counts[rownames(upregulated),], 
                   raw_counts[rownames(downregulated),])
     pheatmap(degs, #T and F are True and False
              cluster_rows = F,
              cluster_cols = F,
              show_rownames = F,
              scale = 'row',
              show_colnames = T)
     #what are the genes that are upregulated
     rownames(upregulated)
     rownames(downregulated)
     
     #exporting the files
     write.csv(upregulated, 'upregulated.csv')
     write.csv(downregulated, 'downregulated.csv')
     write.csv(raw_counts, 'raw_counts.csv')
     
     #Functional Enrichment Analysis
     # Visit https://bioinformatics.sdstate.edu/go/
```
## Results

## Conclusion

## References
