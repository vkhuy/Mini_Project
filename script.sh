#!usr/bin/env bash

# Quality Control for short_reads
mkdir reads_qc
mamba run -n fastp fastp --in1 reads/SRR8278185_1_sub.fastq.gz --in2 reads/SRR8278185_2_sub.fastq.gz --out1 reads_qc/SRR8278185_1_sub.fastq.gz --out2 reads_qc/SRR8278185_2_sub.fastq.gz --unpaired1 reads_qc/SRR8278185_u.fastq.gz --unpaired2 reads_qc/SRR8278185_u.fastq.gz

# Quality Control for long_read
mamba run -n filtlong filtlong --min_length 3000 --keep_percent 90 reads/SRR8278838_sub.fastq.gz | gzip > reads_qc/SRR8278838_sub.fastq.gz
mamba run -n fastqc fastqc reads_qc/SRR8278838_sub.fastq.gz
mv reads_qc/SRR8278838_sub_fastqc.html ./SRR8278838_sub_fastqc.html

# Assemble with trycycler
mamba run -n trycycler trycycler subsample --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir read_subsets

threads=8
mkdir assemblies

mamba run -n trycycler flye --nano-hq read_subsets/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01
mamba run -n trycycler ./miniasm_and_minipolish.sh read_subsets/sample_02.fastq "$threads" > assemblies/assembly_02.gfa && any2fasta assemblies/assembly_02.gfa > assemblies/assembly_02.fasta
mamba run -n trycycler raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa read_subsets/sample_03.fastq > assemblies/assembly_03.fasta

mamba run -n trycycler flye --nano-hq read_subsets/sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04
mamba run -n trycycler ./miniasm_and_minipolish.sh read_subsets/sample_05.fastq "$threads" > assemblies/assembly_05.gfa && any2fasta assemblies/assembly_05.gfa > assemblies/assembly_05.fasta
mamba run -n trycycler raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa read_subsets/sample_06.fastq > assemblies/assembly_06.fasta

mamba run -n trycycler flye --nano-hq read_subsets/sample_07.fastq --threads "$threads" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07
mamba run -n trycycler ./miniasm_and_minipolish.sh read_subsets/sample_08.fastq "$threads" > assemblies/assembly_08.gfa && any2fasta assemblies/assembly_08.gfa > assemblies/assembly_08.fasta
mamba run -n trycycler raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa read_subsets/sample_09.fastq > assemblies/assembly_09.fasta

mamba run -n trycycler flye --nano-hq read_subsets/sample_10.fastq --threads "$threads" --out-dir assembly_10 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa assemblies/assembly_10.gfa && rm -r assembly_10
mamba run -n trycycler ./miniasm_and_minipolish.sh read_subsets/sample_11.fastq "$threads" > assemblies/assembly_11.gfa && any2fasta assemblies/assembly_11.gfa > assemblies/assembly_11.fasta
mamba run -n trycycler raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_12.gfa read_subsets/sample_12.fastq > assemblies/assembly_12.fasta

mkdir assemblies/bad_assembly
mv assemblies/assembly_2.fasta > assemblies/bad_assembly/assembly_2.fasta
mv assemblies/assembly_2.gfa > assemblies/bad_assembly/assembly_2.gfa
mv assemblies/assembly_4.fasta > assemblies/bad_assembly/assembly_4.fasta
mv assemblies/assembly_4.gfa > assemblies/bad_assembly/assembly_4.gfa
mv assemblies/assembly_7.fasta > assemblies/bad_assembly/assembly_7.fasta
mv assemblies/assembly_7.gfa > assemblies/bad_assembly/assembly_7.gfa
mv assemblies/assembly_11.fasta > assemblies/bad_assembly/assembly_11.fasta
mv assemblies/assembly_11.gfa > assemblies/bad_assembly/assembly_11.gfa

mamba run -n trycycler trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir trycycler --distance 0.02

mv trycycler/cluster_003 trycycler/bad_cluster_003

mamba run -n trycycler trycycler reconcile --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dir trycycler/cluster_001 --min_identity 97
mv trycycler/cluster_001/1_contigs/B_Utg1910.fasta trycycler/cluster_001/1_contigs/B_Utg1910.fasta.bad

mamba run -n trycycler trycycler reconcile --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dir trycycler/cluster_002 
mv trycycler/cluster_002/1_contigs/A_contig_3.fasta trycycler/cluster_002/1_contigs/A_contig_3.fasta.bad
mv trycycler/cluster_002/1_contigs/F_Utg1742.fasta trycycler/cluster_002/1_contigs/F_Utg1742.fasta.bad
mv trycycler/cluster_002/1_contigs/G_contig_2.fasta trycycler/cluster_002/1_contigs/G_contig_2.fasta.bad

mamba run -n trycycler trycycler msa --cluster_dir trycycler/cluster_001
mamba run -n trycycler trycycler msa --cluster_dir trycycler/cluster_002

mamba run -n trycycler trycycler partition --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dirs trycycler/cluster_*

# Short_reads polishing with polypolish
mamba run -n polypolish bwa index trycycler.fasta
mamba run -n polypolish bwa mem -t 8 -a trycycler.fasta reads_qc/SRR8278185_1_sub.fastq.gz > alignments_1.sam
mamba run -n polypolish bwa mem -t 8 -a trycycler.fasta reads_qc/SRR8278185_2_sub.fastq.gz > alignments_2.sam
mamba run -n polypolish polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
mamba run -n polypolish polypolish trycycler.fasta filtered_1.sam filtered_2.sam > trycycler_polypolish.fasta
rm *.bwt *.pac *.ann *.amb *.sa *.sam

# Species classification with kraken2
mamba run -n kraken2 kraken2 --db kraken2/database/library/k2_standard_16gb_20230605 trycycler_polypolish.fasta  --report kraken.report --use-names

# Annotation with prokka
mamba run -n prokka prokka --outdir prokka_results --genus helicobacter --species pylori  trycycler_polypolish.fasta

