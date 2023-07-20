#!usr/bin/env bash

# Import environment
bash Modules/Module01.sh


# Quality Control 
mkdir reads_qc
bash Modules/Module02.sh
mkdir qc_report
mv fastp.html qc_report/fastp.html
mv fastp.json qc_report/fastp.json
mv reads_qc/SRR8278838_sub_fastqc.html qc_report/SRR8278838_sub_fastqc.html


# Assemble with trycycler
bash Modules/Module03.sh
# Merge together
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta


# Long_reads polishing with medaka
bash Modules/Module04.sh
# Clean up
mv trycycler/cluster_001/medaka/consensus.fasta trycycler/cluster_001/8_medaka.fasta
rm -r trycycler/cluster_001/medaka trycycler/cluster_001/*.fai trycycler/cluster_001/*.mmi
mv trycycler/cluster_002/medaka/consensus.fasta trycycler/cluster_002/8_medaka.fasta
rm -r trycycler/cluster_002/medaka trycycler/cluster_002/*.fai trycycler/cluster_002/*.mmi 
# Merge together
cat trycycler/cluster_*/8_medaka.fasta > trycycler_medaka.fasta


# Short_reads polishing with polypolish
bash Modules/Module05.sh
rm *.bwt *.pac *.ann *.amb *.sa *.sam

# Move all assembled genomes to assembled folder
mkdir assembled
mv trycycler.fasta assembled/trycycler.fasta
mv trycycler_medaka.fasta assembled/trycycler_medaka.fasta
mv trycycler_medaka_polypolish.fasta assembled/trycycler_medaka_polypolish.fasta

# Species classification with kraken2
bash Modules/Module06.sh


# Annotation with prokka
bash Modules/Module07.sh
mkdir kraken2_report
mv kraken.report kraken2_report/kraken.report
mv SRR8278838.kraken kraken2_report/SRR8278838.kraken
