# #!usr/bin/env bash

# Input
input1=$1 ### Short-read
input2=$2
input3=$3 ### Long-read

# Dowload any2fasta
wget https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta
chmod +x any2fasta 

# Extract the filename without the extension 
filename1=$(basename -- "$input1")
filename1="${filename1%.*.*}"
filename2=$(basename -- "$input2")
filename2="${filename2%.*.*}"
filename3=$(basename -- "$input3")
filename3="${filename3%.*.*}"

# Assemble with trycycler
# Trycycler step 1: subsampling reads
mamba run -n trycycler trycycler subsample --reads reads_qc/${filename3}.fastq.gz --out_dir read_subsets

# Trycycler step 2: assembly
threads=8
mkdir assemblies

for sample in 01 04 07 10; do
    mamba run -n trycycler flye --nano-hq read_subsets/sample_"$sample".fastq --threads "$threads" --out-dir assembly_"$sample" && cp assembly_"$sample"/assembly.fasta assemblies/assembly_"$sample".fasta && cp assembly_"$sample"/assembly_graph.gfa assemblies/assembly_"$sample".gfa && rm -r assembly_"$sample"
done

for sample in 02 05 08 11; do
    mamba run -n trycycler ./miniasm_and_minipolish.sh read_subsets/sample_"$sample".fastq "$threads" > assemblies/assembly_"$sample".gfa && ./any2fasta assemblies/assembly_"$sample".gfa > assemblies/assembly_"$sample".fasta
done

for sample in 03 06 09 12; do
    mamba run -n trycycler raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_"$sample".gfa read_subsets/sample_"$sample".fastq > assemblies/assembly_"$sample".fasta
done
mkdir assemblies/bad_assembly

# Remove bad assemblies
for file in 02 04 07 11; do
    mv assemblies/assembly_"$file".fasta assemblies/bad_assembly/assembly_"$file".fasta
    mv assemblies/assembly_"$file".gfa assemblies/bad_assembly/assembly_"$file".gfa
done

# Trycycler step 3: clustering contigs
mamba run -n trycycler trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir trycycler --distance 0.02

# Visualize clusters tree and found that cluster_003 was bad
java -jar figtree.jar -graphic PNG trycycler/contigs.newick trycycler/clusters.png

# Remove bad clusters
mv trycycler/cluster_003 trycycler/bad_cluster_003

# Remove bad contigs from cluster_001
mv trycycler/cluster_001/1_contigs/B_Utg1910.fasta trycycler/cluster_001/1_contigs/B_Utg1910.fasta.bad

# Remove bad contigs from cluster_002
for file in "A_contig_3" "F_Utg1742" "G_contig_2"; do
    mv trycycler/cluster_002/1_contigs/"$file".fasta trycycler/cluster_002/1_contigs/"$file".fasta.bad
done

# Trycycler step 4: reconciling clusters
mamba run -n trycycler trycycler reconcile --reads reads_qc/${filename3}.fastq.gz --cluster_dir trycycler/cluster_001 --min_identity 97
mamba run -n trycycler trycycler reconcile --reads reads_qc/${filename3}.fastq.gz --cluster_dir trycycler/cluster_002 

# Trycycler step 5: multiple sequence alignment
mamba run -n trycycler trycycler msa --cluster_dir trycycler/cluster_001
mamba run -n trycycler trycycler msa --cluster_dir trycycler/cluster_002

# Trycycler step 6: read partitioning
mamba run -n trycycler trycycler partition --reads reads_qc/${filename3}.fastq.gz --cluster_dirs trycycler/cluster_*

# Trycycler step 7: consensus generation
mamba run -n trycycler trycycler consensus --cluster_dir trycycler/cluster_001
mamba run -n trycycler trycycler consensus --cluster_dir trycycler/cluster_002

echo "Assembly with trycycler: Done"