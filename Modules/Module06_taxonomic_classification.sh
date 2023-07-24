#!usr/bin/env bash

# Taxonomic classification with kraken2
# Dowload and unzip database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20230605.tar.gz
tar -xf k2_standard_16gb_20230605.tar.gz -C kraken2/

mamba run -n kraken2 kraken2 --db kraken2/k2_standard_16gb_20230605 trycycler_medaka_polypolish.fasta  --report kraken.report --use-names > SRR8278838.kraken
