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
<insert code> 
```
## Results

## Conclusion

## References
