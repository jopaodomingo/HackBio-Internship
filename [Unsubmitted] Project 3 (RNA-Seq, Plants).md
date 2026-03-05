# Project 3 (RNA-Seq Plants)
## Introduction
In this project, analysis of the <i>Arabidopsis thaliana</i> leaf's response to UV stress in the vasculature was conducted using Bash and RStudio. The objectives include the following: 
1) To perform differential expression analysis in the vasculature of UV treated <i>Arabidopsis thaliana</i> tissues with Water (control) and UV-C light
2) To express the results via volcano plot and heat map using RStudio
3) To generate a list of upregulated and downregulated genes as a .csv file

## Methodology
### Dataset
The data set consisted of 6 samples that were obtained from SRA-Explorer. These samples consisted of 3 replicates that were exposed to water (SRR12808527, SRR12808528, SRR12808529) and UV-C light (SRR12808497, SRR12808498, SRR12808499). A series of pre-processing steps were conducted as follows:

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
The data is then saved in a txt file named "counts.txt" and uploaded to RStudio, where Differential Sequencing Analysis (DSeq2) is applied. 
### Differential Expression Analysis
```
#Activation of library
library(DESeq2) 

#Input of data file "counts.txt" and "metadata.txt" containing the sample names and whether they are control samples or exposed to UV-C
c_e_count <- read.delim('counts.txt', header = T)
c_e_meta <- read.delim('metadata.txt', header = T, stringsAsFactors = T) 

#Consider only columns SRR12808497 to SRR12808529
raw_counts <- c_e_count[, c("SRR12808497", "SRR12808498", "SRR12808499", "SRR12808527", "SRR12808528", "SRR12808529")]

#Row numbers are converted to the gene IDs
rownames(raw_counts) <- c_e_count$Geneid

#Create desqdataset; countData is the actual data; colData is the metadata for differential expression analysis 
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = c_e_meta,
                              design = ~stress_exposure)

#Perform the differential expression analysis
dds <- DESeq(dds)

#Assign the final results to variable final_res
final_res <- results(dds)

#Check distribution of p-values
plot(density(x = na.omit(final_res$pvalue)))

#Create Volcano plot of differentially expressed genes
plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$padj),  #x is the magnitude of change, y is the distribution of p values 
     cex = 0.25,#size of characters
     pch = 19, #shape used 
     col = 'grey', #color
     ylim = c(0,20), #y-range 
     ylab = 'Adjusted P-Value', #label
     xlab = 'Log2 FC') #label

     abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)
     
     #Highlight upregulated genes in UV-C compared to control samples 
     upregulated <- subset(final_res, padj < 0.05 & log2FoldChange > 2) #log2FoldChange is change in magnitude
     points(upregulated$log2FoldChange,
            y = -log10(upregulated$padj), 
            cex = 0.35,
            pch = 19,
            col = 'salmon')
     
     #Highlight downregulated genes in UV-C compared to control samples
     downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -2)
     points(downregulated$log2FoldChange,
            y = -log10(downregulated$padj), 
            cex = 0.35,
            pch = 19,
            col = 'lightblue')
     #Insert Title
     mtext('Volcano plot of gene regulation in UV-C stressed Arabidopsis thaliana')
```
A heat map of the upregulated and downregulated genes is then generated as follows:
### Heat Map Generation
```
#Activate library
library(pheatmap)

#Set degs as input data
degs <- rbind(raw_counts[rownames(upregulated),], 
        raw_counts[rownames(downregulated),])
#Heat map
pheatmap(degs, #T and F are True and False
        cluster_rows = F,
        cluster_cols = F,
        show_rownames = F,
        scale = 'row',
        show_colnames = T)
```
The list of upregulated and downregulated genes are then exported as .csv files:
### Exportation of list of upregulated and downregulated genes
```
write.csv(upregulated, 'upregulated.csv')
write.csv(downregulated, 'downregulated.csv')
```
## Results
The distribution of p-values from the results of the differential expression analysis is shown below: 
<img width="725" height="585" alt="P" src="https://github.com/user-attachments/assets/4c05c74b-79bf-48df-b3db-49a3da71b0eb" />

The generated volcano plot of the differential expression analysis is shown below. Note that the statistically significant (p < 0.05) upregulated genes in UV-C exposed samples are highlighted as pink and the statistically significant downregulated genes in UV-C exposed samples are highlighted as blue. Note that log2FC refers to the log2 Fold Change or the magnitude of change in gene expression. 
<img width="725" height="585" alt="Volcano Plot" src="https://github.com/user-attachments/assets/6c7b4186-0551-4c3d-9fb7-5621853ac3cd" />

The produced heat map for each sample ID is presented below. Recall that samples exposed to UV-C light include SRR12808497, SRR12808498, and SRR12808499, while control samples exposed to water include SRR12808527, SRR12808528, and SRR12808529. It can be noticed that a larger proportion of genes are upregulated in UV-C exposed samples compared to their downregulated counterparts. 
<img width="725" height="585" alt="Heat_Map" src="https://github.com/user-attachments/assets/a949538e-1acf-4e03-aa1c-eae0bedc713f" />

## Recommendations 
Further studies of this project may consider Functional Enrichment Analysis to determine the top 100 upregulated genes in the UV-C exposed samples. 
