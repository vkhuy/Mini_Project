#!usr/bin/env bash

### Short_reads polishing with polypolish

# BWA the assembly with short reads
mamba run -n polypolish bwa index trycycler_medaka.fasta
mamba run -n polypolish bwa mem -t 8 -a trycycler_medaka.fasta reads_qc/SRR8278185_1_sub.fastq.gz > alignments_1.sam
mamba run -n polypolish bwa mem -t 8 -a trycycler_medaka.fasta reads_qc/SRR8278185_2_sub.fastq.gz > alignments_2.sam

# Polishing
mamba run -n polypolish polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
mamba run -n polypolish polypolish trycycler_medaka.fasta filtered_1.sam filtered_2.sam > trycycler_medaka_polypolish.fasta
