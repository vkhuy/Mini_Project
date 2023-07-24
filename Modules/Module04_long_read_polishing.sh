#!usr/bin/env bash

# Long_reads polishing with medaka in 2 clusters
mamba run -n medaka medaka_consensus -i trycycler/cluster_001/4_reads.fastq -d trycycler/cluster_001/7_final_consensus.fasta -o trycycler/cluster_001/medaka

mamba run -n medaka medaka_consensus -i trycycler/cluster_002/4_reads.fastq -d trycycler/cluster_002/7_final_consensus.fasta -o trycycler/cluster_002/medaka

