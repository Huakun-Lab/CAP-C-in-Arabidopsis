### PLEASE PUT YOUR CAPC DATA INTO "01-CAPC" DIRECTORY AND RENAME TO "CAPC_R1(R2).fq.gz"!! ###

PIGZ=/PATH/TO/YOUR/PIGZ/pigz
SEQPREP=/PATH/TO/YOUR/SEQPREP/SeqPrep
CUTADAPT=/PATH/TO/YOUR/CUTADAPT/cutadapt
PYTHON=/PATH/TO/YOUR/PYTHON/python

FORWARD_ADAPTER="ACGCGATATCTTATCTGACT" 
REVERSE_ADAPTER="AGTCAGATAAGATATCGCGT"

mkdir -p 20-split
mkdir -p 21-merge
mkdir -p 02-cutadapt
mkdir -p 03-after_filter
mkdir -p 04-fastq/CAPC

## Step 1 Split raw fastq file ##

for FILE in `ls 01-CAPC`
do
	NAME=`echo $FILE | sed 's/\.fq.gz//'`
	zcat 01-CAPC/$FILE | split -a 3 -d -l 20000000 - 20-split/$NAME\_ 
done

## Step 2 Zip Split fastq file ##

for FILE in `ls 20-split`
do
	mv 20-split/$FILE 20-split/$FILE\.fastq
	$PIGZ -p 60 20-split/$FILE\.fastq 
done

## Step 3 Merge R1 and R2 ##

for i in `ls 20-split/CAPC_R1* | sed 's/_/\t/g' | sed 's/\./\t/g' | awk '{print $3}'`
do
	$SEQPREP -f 20-split/CAPC_R1_$i\.fastq.gz \
	-r 20-split/CAPC_R2_$i\.fastq.gz \
	-s 21-merge/CAPC_merged_$i\.fastq.gz \
	-1 21-merge/CAPC_unmerged_R1_$i\.fastq.gz \
	-2 21-merge/CAPC_unmerged_R2_$i\.fastq.gz
done

## Step 4 Remove Linker ##

for i in `ls 20-split/$REP\_R1* | sed 's/_/\t/g' | sed 's/\./\t/g' | awk '{print $3}'`
do
	$CUTADAPT -n 1 --overlap 10 -a forward=$FORWARD_ADAPTER -a reverse=$REVERSE_ADAPTER -o 02-cutadapt/CAPC_merged_$i\_T1.fastq.gz 21-merge/CAPC_merged_$i\.fastq.gz 
	$CUTADAPT -n 1 --overlap 10 -g forward=$FORWARD_ADAPTER -g reverse=$REVERSE_ADAPTER -o 02-cutadapt/CAPC_merged_$i\_T2.fastq.gz 21-merge/CAPC_merged_$i\.fastq.gz 
	$CUTADAPT -n 1 --overlap 10 -a forward=$FORWARD_ADAPTER -a reverse=$REVERSE_ADAPTER -o 02-cutadapt/CAPC_unmerged_$i\_A1.fastq.gz --mask-adapter 21-merge/CAPC_unmerged_R1_$i\.fastq.gz 
	$CUTADAPT -n 1 --overlap 10 -g forward=$FORWARD_ADAPTER -g reverse=$REVERSE_ADAPTER -o 02-cutadapt/CAPC_unmerged_$i\_G1.fastq.gz --mask-adapter 21-merge/CAPC_unmerged_R1_$i\.fastq.gz 
	$CUTADAPT -n 1 --overlap 10 -a forward=$FORWARD_ADAPTER -a reverse=$REVERSE_ADAPTER -o 02-cutadapt/CAPC_unmerged_$i\_A2.fastq.gz --mask-adapter 21-merge/CAPC_unmerged_R2_$i\.fastq.gz 
	$CUTADAPT -n 1 --overlap 10 -g forward=$FORWARD_ADAPTER -g reverse=$REVERSE_ADAPTER -o 02-cutadapt/CAPC_unmerged_$i\_G2.fastq.gz --mask-adapter 21-merge/CAPC_unmerged_R2_$i\.fastq.gz 
done

## Step 5 Sum Split and filted ##

for i in `ls 20-split/CAPC_R1* | sed 's/_/\t/g' | sed 's/\./\t/g' | awk '{print $3}'`
do
	$PYTHON src/Rewrite_FilterN.py -m 20 -n 10 \
	-m1 02-cutadapt/CAPC_merged_$i\_T1.fastq.gz \
	-m2 02-cutadapt/CAPC_merged_$i\_T2.fastq.gz \
	-A1 02-cutadapt/CAPC_unmerged_$i\_A1.fastq.gz \
	-G1 02-cutadapt/CAPC_unmerged_$i\_G1.fastq.gz \
	-A2 02-cutadapt/CAPC_unmerged_$i\_A2.fastq.gz \
	-G2 02-cutadapt/CAPC_unmerged_$i\_G2.fastq.gz \
	-prefix 03-after_filter/CAPC_merged_$i & 
done

/data/software/01-anaconda3/envs/huakun_py2/bin/pigz -p 60 03-after_filter/*

## Step 6 Combine all fastq ##

cat 03-after_filter/CAPC*{OK,Short}_R1*.fastq.gz > 04-fastq/CAPC/CAPC_merged_R1.fastq.gz
cat 03-after_filter/CAPC*{OK,Short}_R2*.fastq.gz > 04-fastq/CAPC/CAPC_merged_R2.fastq.gz

