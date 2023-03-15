# CAP-C-in-Arabidopsis
The following instructions describe the pipeline used to process Arabidopsis CAP-C libraries in the manuscript. 

## Dependencies and Prerequisites
The following software are required to run the Bash pipeline

- pigz 
- SeqPrep
- cutadapt
- python 3.7
- HiC-Pro (>=3.0.0) 


## File storage structure
- 01-CAPC
	- CAPC\_R1.fq.gz
	- CAPC\_R2.fq.gz
- 02-cutadapt
	- CAPC\_merged\_001\_T1.fastq.gz
	- CAPC\_merged\_001\_T2.fastq.gz
	- CAPC\_unmerged\_001\_A1.fastq.gz
	- CAPC\_unmerged\_001\_G1.fastq.gz
	- CAPC\_unmerged\_001\_A2.fastq.gz
	- CAPC\_unmerged\_001\_G2.fastq.gz
	- ......
- 03-after_filter
	- CAPC\*{OK,Short}\_R1\*.fastq.gz
	- CAPC\*{OK,Short}\_R2\*.fastq.gz
- 04-fastq
	- CAPC\_merged\_R1.fastq.gz
	- CAPC\_merged\_R2.fastq.gz
- 20-split
	- CAPC\_R1\_001.fastq.gz
	- CAPC\_R1\_002.fastq.gz
	- ......
- 21-merge
	- CAPC\_merged\_001\.fastq.gz
	- CAPC\_unmerged\_R1\_001\.fastq.gz
	- CAPC\_unmerged\_R2\_001\.fastq.gz
	- ......

When the Bash Script finished. You can remove folder `02-cutadapt`, `03-after_filter`, `20-split`, `21-merge`
The clean fastq was storage in `04-fastq`.

## General Instructions to use the bash script