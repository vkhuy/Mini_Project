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

### Setup environment:
# Import environment
bash Modules/setup_environment.sh

### Module01:
# Quality Control 
mkdir reads_qc
bash Modules/Module01_quality_control.sh $input1 $input2 $input3

### Module02:
# Assembly with trycycler
bash Modules/Module02_assembly.sh $input1 $input2 $input3
# Merge together
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta

### Module03:
# Long_reads polishing with medaka
bash Modules/Module03_long_read_polishing.sh
# Clean up
mv trycycler/cluster_001/medaka/consensus.fasta trycycler/cluster_001/8_medaka.fasta
rm -r trycycler/cluster_001/medaka trycycler/cluster_001/*.fai trycycler/cluster_001/*.mmi
mv trycycler/cluster_002/medaka/consensus.fasta trycycler/cluster_002/8_medaka.fasta
rm -r trycycler/cluster_002/medaka trycycler/cluster_002/*.fai trycycler/cluster_002/*.mmi 
# Merge together
cat trycycler/cluster_*/8_medaka.fasta > trycycler_medaka.fasta

### Module04:
# Short_reads polishing with polypolish
bash Modules/Module04_short_read_polishing.sh $input1 $input2 
rm *.bwt *.pac *.ann *.amb *.sa *.sam

# Move all assembled genomes to assembled folder
mkdir assembled
mv trycycler.fasta assembled/trycycler.fasta
mv trycycler_medaka.fasta assembled/trycycler_medaka.fasta
mv trycycler_medaka_polypolish.fasta assembled/trycycler_medaka_polypolish.fasta
echo "Assembly: Done"

### Module05:
Taxonomic classification with kraken2
mkdir kraken2
bash Modules/Module05_taxonomic_classification.sh $input3

## Module06:
# Annotation with prokka
bash Modules/Module06_Annotation.sh


echo "Pipeline is finished"