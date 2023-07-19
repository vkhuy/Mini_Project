#!usr/bin/env bash

# Quality Control for short_reads: remove apdapters and trim off low-quality bases
mamba run -n fastp fastp --in1 reads/SRR8278185_1_sub.fastq.gz --in2 reads/SRR8278185_2_sub.fastq.gz --out1 reads_qc/SRR8278185_1_sub.fastq.gz --out2 reads_qc/SRR8278185_2_sub.fastq.gz --unpaired1 reads_qc/SRR8278185_u.fastq.gz --unpaired2 reads_qc/SRR8278185_u.fastq.gz

# Quality Control for long_read: remove short reads to improve N50 and throw out the worst 10% of reads
mamba run -n filtlong filtlong --min_length 3000 --keep_percent 90 reads/SRR8278838_sub.fastq.gz | gzip > reads_qc/SRR8278838_sub.fastq.gz
mamba run -n fastqc fastqc reads_qc/SRR8278838_sub.fastq.gz
