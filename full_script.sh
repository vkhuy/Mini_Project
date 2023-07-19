#!usr/bin/env bash

# Import environment
mamba env create --file envs/*.yml


# Quality Control for short_reads
mkdir reads_qc
mamba run -n fastp fastp --in1 reads/SRR8278185_1_sub.fastq.gz --in2 reads/SRR8278185_2_sub.fastq.gz --out1 reads_qc/SRR8278185_1_sub.fastq.gz --out2 reads_qc/SRR8278185_2_sub.fastq.gz --unpaired1 reads_qc/SRR8278185_u.fastq.gz --unpaired2 reads_qc/SRR8278185_u.fastq.gz

# Quality Control for long_read
mamba run -n filtlong filtlong --min_length 3000 --keep_percent 90 reads/SRR8278838_sub.fastq.gz | gzip > reads_qc/SRR8278838_sub.fastq.gz
mamba run -n fastqc fastqc reads_qc/SRR8278838_sub.fastq.gz
mv reads_qc/SRR8278838_sub_fastqc.html ./SRR8278838_sub_fastqc.html


# Assemble with trycycler
# Trycycler step 1: subsampling reads
mamba run -n trycycler trycycler subsample --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir read_subsets

# Trycycler step 2: assembly
threads=8
mkdir assemblies

for sample in 01 04 07 10; do
    mamba run -n trycycler flye --nano-hq read_subsets/sample_"$sample".fastq --threads "$threads" --out-dir assembly_"$sample" && cp assembly_"$sample"/assembly.fasta assemblies/assembly_"$sample".fasta && cp assembly_"$sample"/assembly_graph.gfa assemblies/assembly_"$sample".gfa && rm -r assembly_"$sample"
done

for sample in 02 05 08 11; do
    mamba run -n trycycler ./miniasm_and_minipolish.sh read_subsets/sample_"$sample".fastq "$threads" > assemblies/assembly_"$sample".gfa && any2fasta assemblies/assembly_"$sample".gfa > assemblies/assembly_"$sample".fasta
done

for sample in 03 06 09 12; do
    mamba run -n trycycler raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_"$sample".gfa read_subsets/sample_"$sample".fastq > assemblies/assembly_"$sample".fasta
done
mkdir assemblies/bad_assembly
for file in 2 4 7 11; do
    mv assemblies/assembly_"$file".fasta assemblies/bad_assembly/assembly_"$file".fasta
    mv assemblies/assembly_"$file".gfa assemblies/bad_assembly/assembly_"$file".gfa
done

# Trycycler step 3: clustering contigs
mamba run -n trycycler trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir trycycler --distance 0.02

mv trycycler/cluster_003 trycycler/bad_cluster_003

# Trycycler step 4: reconciling clusters
mamba run -n trycycler trycycler reconcile --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dir trycycler/cluster_001 --min_identity 97
mv trycycler/cluster_001/1_contigs/B_Utg1910.fasta trycycler/cluster_001/1_contigs/B_Utg1910.fasta.bad

mamba run -n trycycler trycycler reconcile --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dir trycycler/cluster_002 

for file in "A_contig_3" "F_Utg1742" "G_contig_2"; do
    mv trycycler/cluster_002/1_contigs/"$file".fasta trycycler/cluster_002/1_contigs/"$file".fasta.bad
done

# Trycycler step 5: multiple sequence alignment
mamba run -n trycycler trycycler msa --cluster_dir trycycler/cluster_001
mamba run -n trycycler trycycler msa --cluster_dir trycycler/cluster_002

# Trycycler step 6: read partitioning
mamba run -n trycycler trycycler partition --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dirs trycycler/cluster_*

# Trycycler step 7: consensus generation
mamba run -n trycycler trycycler consensus --cluster_dir trycycler/cluster_001
mamba run -n trycycler trycycler consensus --cluster_dir trycycler/cluster_002

# Merge together
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler.fasta


# Long_reads polishing with medaka
mamba run -n medaka medaka_consensus -i trycycler/cluster_001/4_reads.fastq -d trycycler/cluster_001/7_final_consensus.fasta -o trycycler/cluster_001/medaka
mv trycycler/cluster_001/medaka/consensus.fasta trycycler/cluster_001/8_medaka.fasta
rm -r trycycler/cluster_001/medaka trycycler/cluster_001/*.fai trycycler/cluster_001/*.mmi

mamba run -n medaka medaka_consensus -i trycycler/cluster_002/4_reads.fastq -d trycycler/cluster_002/7_final_consensus.fasta -o trycycler/cluster_002/medaka
mv trycycler/cluster_002/medaka/consensus.fasta trycycler/cluster_002/8_medaka.fasta
rm -r trycycler/cluster_002/medaka trycycler/cluster_002/*.fai trycycler/cluster_002/*.mmi 

cat trycycler/cluster_*/8_medaka.fasta > trycycler_medaka.fasta


# Short_reads polishing with polypolish
mamba run -n polypolish bwa index trycycler_medaka.fasta
mamba run -n polypolish bwa mem -t 8 -a trycycler_medaka.fasta reads_qc/SRR8278185_1_sub.fastq.gz > alignments_1.sam
mamba run -n polypolish bwa mem -t 8 -a trycycler_medaka.fasta reads_qc/SRR8278185_2_sub.fastq.gz > alignments_2.sam
mamba run -n polypolish polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
mamba run -n polypolish polypolish trycycler_medaka.fasta filtered_1.sam filtered_2.sam > trycycler_medaka_polypolish.fasta
rm *.bwt *.pac *.ann *.amb *.sa *.sam


# Species classification with kraken2
mamba run -n kraken2 kraken2 --db kraken2/database/library/k2_standard_16gb_20230605 trycycler_medaka_polypolish.fasta  --report kraken.report --use-names


# Annotation with prokka
mamba run -n prokka prokka --outdir prokka_results --genus helicobacter --species pylori  trycycler_medaka_polypolish.fasta

