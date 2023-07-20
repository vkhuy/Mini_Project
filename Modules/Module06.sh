#!usr/bin/env bash

# Species classification with kraken2
mamba run -n kraken2 kraken2 --db kraken2/database/library/k2_standard_16gb_20230605 trycycler_medaka_polypolish.fasta  --report kraken.report --use-names > SRR8278838.kraken
