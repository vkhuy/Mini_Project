#!usr/bin/env bash

# Import environment
bash Module01.sh


# Quality Control 
mkdir reads_qc
bash Module02.sh
mv reads_qc/SRR8278838_sub_fastqc.html ./SRR8278838_sub_fastqc.html


# Assemble with trycycler
bash Module03.sh
# Merge together
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta


# Long_reads polishing with medaka
bash Module04.sh
# Clean up
mv trycycler/cluster_001/medaka/consensus.fasta trycycler/cluster_001/8_medaka.fasta
rm -r trycycler/cluster_001/medaka trycycler/cluster_001/*.fai trycycler/cluster_001/*.mmi
mv trycycler/cluster_002/medaka/consensus.fasta trycycler/cluster_002/8_medaka.fasta
rm -r trycycler/cluster_002/medaka trycycler/cluster_002/*.fai trycycler/cluster_002/*.mmi 
# Merge together
cat trycycler/cluster_*/8_medaka.fasta > trycycler_medaka.fasta


# Short_reads polishing with polypolish
bash Module05.sh
rm *.bwt *.pac *.ann *.amb *.sa *.sam


# Species classification with kraken2
bash Module06.sh


# Annotation with prokka
bash Module07.sh
