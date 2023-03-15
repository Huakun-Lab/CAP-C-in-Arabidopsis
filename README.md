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
	- CAPC/CAPC\_merged\_R1.fastq.gz
	- CAPC/CAPC\_merged\_R2.fastq.gz
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
This bash script is used to clean the raw sequencing fastq, then you can use the clean fastq to mapping.
### Step1 Split raw fastq file
Since the software SeqPrep used in the third step generally only performs single-threaded calculations, the original FASTQ files should be split into smaller FASTQ files first.
### Step 2 Zip Split fastq file
In order to reduce storage usage, compress the split FASTQ files. We used pigz, which can perform multi-threaded compression on sequencing files, speeding up data processing.
### Step 3 Merge R1 and R2 
Use SeqPrep to merge the FASTQ reads.
### Step 4 Remove Linker
Remove the linker sequences from the reads and regenerate the FASTQ files.
### Step 5 Sum Split and filted
Filter the regenerated FASTQ files.
### Step 6 Combine all fastq
Reassemble the filtered FASTQ files into clean FASTQ files for the next step of mapping.

## After the bash script
Now we have fastq data without linkers. Then we use `HiC-Pro` to mapping. We provide our HiC-Pro config file `CAPC_Arabidopsis_config-hicpro.txt`, and use the command `hicpro -i 04-fastq -o 05-rmdup_rmmulti_Seperate_HiCPro_result -c Arabidopsis_DNase_config-hicpro.txt `.