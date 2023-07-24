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

### Module01:
# Import environment
bash Modules/Module01_import_environment.sh

### Module02:
# Quality Control 
mkdir reads_qc
bash Modules/Module02_quality_control.sh $input1 $input2 $input3

### Module03:
# Assembly with trycycler
bash Modules/Module03_assembly.sh $input1 $input2 $input3
# Merge together
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta

### Module04:
# Long_reads polishing with medaka
bash Modules/Module04_long_read_polishing.sh
# Clean up
mv trycycler/cluster_001/medaka/consensus.fasta trycycler/cluster_001/8_medaka.fasta
rm -r trycycler/cluster_001/medaka trycycler/cluster_001/*.fai trycycler/cluster_001/*.mmi
mv trycycler/cluster_002/medaka/consensus.fasta trycycler/cluster_002/8_medaka.fasta
rm -r trycycler/cluster_002/medaka trycycler/cluster_002/*.fai trycycler/cluster_002/*.mmi 
# Merge together
cat trycycler/cluster_*/8_medaka.fasta > trycycler_medaka.fasta

### Module05:
# Short_reads polishing with polypolish
bash Modules/Module05_short_read_polishing.sh $input1 $input2 
rm *.bwt *.pac *.ann *.amb *.sa *.sam

# Move all assembled genomes to assembled folder
mkdir assembled
mv trycycler.fasta assembled/trycycler.fasta
mv trycycler_medaka.fasta assembled/trycycler_medaka.fasta
mv trycycler_medaka_polypolish.fasta assembled/trycycler_medaka_polypolish.fasta
echo "Assembly: Done"

### Module06:
Taxonomic classification with kraken2
mkdir kraken2
bash Modules/Module06_taxonomic_classification.sh $input3

## Module07:
# Annotation with prokka
bash Modules/Module07_Annotation.sh


echo "Pipeline is finished"