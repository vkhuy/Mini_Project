#!usr/bin/env bash

# Input
input1=$1 ### Short-read
input2=$2
input3=$3 ### Long-read

# Extract the filename without the extension 
filename1=$(basename -- "$input1")
filename1="${filename1%.*.*}"
filename2=$(basename -- "$input2")
filename2="${filename2%.*.*}"
filename3=$(basename -- "$input3")
filename3="${filename3%.*.*}"

# Quality Control for short_reads: remove apdapters and trim off low-quality bases
mamba run -n fastp fastp --in1 $input1 --in2 $input2 --out1 reads_qc/${filename1}.fastq.gz --out2 reads_qc/${filename2}.fastq.gz --unpaired1 reads_qc/unpair.fastq.gz --unpaired2 reads_qc/unpair.fastq.gz

# Quality Control for long_read: remove short reads to improve N50 and throw out the worst 10% of reads
mamba run -n filtlong filtlong --min_length 3000 --keep_percent 90 $input3 | gzip > reads_qc/${filename3}.fastq.gz
mamba run -n fastqc fastqc reads_qc/${filename3}.fastq.gz

mkdir qc_report
mv fastp.html qc_report/fastp.html
mv fastp.json qc_report/fastp.json
mv reads_qc/${filename3}_fastqc.html qc_report/${filename3}_fastqc.html
mv reads_qc/${filename3}_fastqc.zip qc_report/${filename3}_fastqc.zip

echo "Quality Control: Done"
