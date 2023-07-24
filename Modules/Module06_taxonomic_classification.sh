#!usr/bin/env bash

# Input
input3=$1 ### Long-read

# Extract the filename without the extension 
filename3=$(basename -- "$input3")
filename3="${filename3%.*.*}"

### Taxonomic classification with kraken2
# Dowload and unzip database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20230605.tar.gz
tar -xf k2_standard_16gb_20230605.tar.gz -C kraken2/

mamba run -n kraken2 kraken2 --db kraken2/k2_standard_16gb_20230605 assembled/trycycler_medaka_polypolish.fasta  --report kraken.report --use-names > ${filename3}.kraken

# Move reports
mkdir kraken2_report
mv kraken.report kraken2_report/kraken.report
mv ${filename3}.kraken kraken2_report/${filename3%.*.*}.kraken

echo "Taxonomic classification: Done"