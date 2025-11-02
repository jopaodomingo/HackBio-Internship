# Project 2 - WGS Microbes
In this project, we will be using Whole-Genome Sequencing (WGS) to analyze data from 50 Bacterial isolates collected during the 2017-2018 South African Outbreak. Specifically, the objectives of this project are as follows:
1) Confirm the identity of the organism
2) Determine the antimicrobial resistance (AMR) profile of these pathogens.
3) Detect if there might be a toxin that is accelerating the death rate
4) Suggest antibiotics/treatment options that could be used to imaange the cases

## Dataset
The data set was obtained from a subset of 100 samples. Downsizing was done due to time constraints. Data was downloading from the following site: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR270/016/SRR27013316/ using the curl command and stored in a directory entitled "raw_reads". Samples used include SRR27013147, SRR27013245, SRR27013246,..., SRR27013292, SRR27013293. 

## FastQC implementation
> nano run_fastQC.sh
```
#!/bin/bash 

#Creation of directories 
mkdir qc_report 

#For_loop

for filename in raw_reads/*.fastq.gz; do 
	fastqc "$filename" -o qc_report/
done
```
Issues in the per base content of the samples (varying %base content in bp positions 1-15) suggested the need for trimming using FastP:
## FastP implementation
> nano run_fastp.sh
```
#!/bin/bash 
#Creation of directory
mkdir trimmed_reads

#For loop for fastp
for i in {147} {245..293}; do
gene_1="SRR27013${i}_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_1"
gene_2="SRR27013${i}_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_2"
gene_name="SRR27013${i}_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017"

fastp \
    	-i "raw_reads/${gene_1}.fastq.gz" \
    	-I "raw_reads/${gene_2}.fastq.gz" \
    	-o "trimmed_reads/${gene_1}_trimmed.fastq.gz" \
    	-O "trimmed_reads/${gene_2}_trimmed.fastq.gz" \
    	--html "trimmed_reads/${gene_name}_fastp.html" #generate an HTML report for each sample 
done
```
After trimming, de novo genome assembly was done using spades.py. 
## Assembly
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
Finally, to study antimcriobial resistance, AMR profile was obtained using abricate:
## AMR Profile
> nano run_AMR.sh
```
#!/bin/bash

#For loop for abricate
for i in {147} {245..293}; do
	mkdir "assembly/SRR27013${i}/AMR"
  abricate "assembly/SRR27013${i}/contigs.fasta" \
> "assembly/SRR27013${i}/AMR/amr_tab.tab"
done

```
## Results
<img width="1895" height="958" alt="image" src="https://github.com/user-attachments/assets/2663f764-6c09-41df-afdb-847597a1f2a4" />
Using BLAST interface, the genome involved was determined to be Listeria monocytogenes. 

