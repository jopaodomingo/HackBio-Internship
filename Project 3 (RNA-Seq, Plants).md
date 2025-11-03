# Project 3 (RNA-Seq Plants)
## Introduction
In this project, analysis of the <i>Arabidopsis thaliana</i> leaf's response to UV stress in the vasculature was conducted using Bash and RStudio. The objectives include the following: 
1) To perform differential expression analysis in the vasculature of UV treated <i>Arabidopsis thaliana</i> tissues with Water (control) and UV-C light
2) To determine the top 100 differentially expressed genes of the two groups of <i>Arabidopsis thaliana</i> tissues
3) To determine the top 5 enriched pathways of the two groups of <i>Arabidopsis thaliana</i> tissues

## Methodology
### Dataset
The data set consisted of 6 samples that were obtained from SRA-Explorer. These samples consisted of 3 replicates that were exposed to water (SRR12808527, SRR12808528, SRR12808529) and UV-C light (SRR12808497, SRR12808498, SRR12808499).

A series of pre-processing steps were conducted as follows 
### FastQC implementation
```
#!/bin/bash 
mkdir qc

for filename in clean_files/*.gz; do 
  fastqc -o qc/ $filename
  gunzip
done
```
<img width="1506" height="813" alt="image" src="https://github.com/user-attachments/assets/9dd9de2b-0e41-40a0-8f12-6f84fc10f6d1" />
Issues in the per base sequence content, sequence duplication levels, and overrepresented sequences of the samples suggested the need for trimming via FastP.
### FastP implementation
> nano run_fastp.sh
```
#!/bin/bash 
#Creation of directory
mkdir trim

#For loop
for filename in clean_files/*.fastq; do

```
After trimming, reference genome mapping was done using STAR
### Reference genome mapping using STAR
> nano run_assembly.sh
```
#!/bin/bash
#Creation of directory
mkdir assembly

#For loop for assembly
for i in {147} {245...293}; do
	gene_1="SRR27013${i}_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_1"
	gene_2="SRR27013${i}_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_2"
	mkdir "assembly/SRR27013${i}"
	spades.py \
		--phred-offset 33 \
    		-1 "trimmed_reads/${gene_1}_trimmed.fastq.gz" \
    		-2 "trimmed_reads/${gene_2}_trimmed.fastq.gz" \
    		-o "assembly/SRR27013${i}"
done

```
The Bandage software was used to open the respectived .fastg files of each assembled genome. A random node with less than 200 bp was selected and inputted into the BLAST interface (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=MegaBlast&PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=) to determine any matching reference genomes to identify the genome species.
<center><figure>
	<img width="1920" height="1032" alt="image" src="https://github.com/user-attachments/assets/9de5ffb8-952b-471e-b515-df1ff10c5004" />
	<figcaption><center><i>Assembled genome for SRR27013260 using Bandage software</i></center></figcaption>
</figure></center>
<br> <br>
<center><figure>
	<img width="1895" height="958" alt="image" src="https://github.com/user-attachments/assets/2663f764-6c09-41df-afdb-847597a1f2a4" />
	<figcaption><center><i>Blast output for SRR27013147</i></center></figcaption>
</figure></center>
<br> <br>

Finally, to study antimcriobial resistance, AMR profile was obtained using abricate:
### AMR Profile
> nano run_AMR.sh
```
#!/bin/bash

#For loop for abricate
for i in {147} {245..293}; do
	mkdir "assembly/SRR27013${i}/AMR"
  abricate "assembly/SRR27013${i}/contigs.fasta" \
> "assembly/SRR27013${i}/AMR/amr_tab.tab"
done
## Results

## Conclusion

## References
