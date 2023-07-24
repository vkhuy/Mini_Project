#!usr/bin/env bash

# Input
input1=$1 ### Short-read
input2=$2

# Extract the filename without the extension 
filename1=$(basename -- "$input1")
filename1="${filename1%.*.*}"
filename2=$(basename -- "$input2")
filename2="${filename2%.*.*}"

### Short_reads polishing with polypolish

# BWA the assembly with short reads
mamba run -n polypolish bwa index trycycler_medaka.fasta
mamba run -n polypolish bwa mem -t 8 -a trycycler_medaka.fasta reads_qc/${filename1}.fastq.gz > alignments_1.sam
mamba run -n polypolish bwa mem -t 8 -a trycycler_medaka.fasta reads_qc/${filename2}.fastq.gz > alignments_2.sam

# Polishing
mamba run -n polypolish polypolish trycycler_medaka.fasta alignments_1.sam alignments_2.sam > trycycler_medaka_polypolish.fasta

echo "Short_reads polishing: Done"


