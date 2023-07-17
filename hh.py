

mkdir reads_qc
fastp --in1 reads/SRR8278185_1_sub.fastq.gz --in2 reads/SRR8278185_2_sub.fastq.gz --out1 reads_qc/SRR8278185_1_sub.fastq.gz --out2 reads_qc/SRR8278185_2_sub.fastq.gz --unpaired1 reads_qc/SRR8278185_u.fastq.gz --unpaired2 reads_qc/SRR8278185_u.fastq.gz


filtlong --min_length 3000 --keep_percent 90 reads/SRR8278838_sub.fastq.gz | gzip > reads_qc/SRR8278838_sub.fastq.gz

trycycler subsample --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir read_subsets

threads=8
mkdir assemblies

flye --nano-hq read_subsets/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01
./miniasm_and_minipolish.sh read_subsets/sample_02.fastq "$threads" > assemblies/assembly_02.gfa && any2fasta assemblies/assembly_02.gfa > assemblies/assembly_02.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa read_subsets/sample_03.fastq > assemblies/assembly_03.fasta

flye --nano-hq read_subsets/sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04
./miniasm_and_minipolish.sh read_subsets/sample_05.fastq "$threads" > assemblies/assembly_05.gfa && any2fasta assemblies/assembly_05.gfa > assemblies/assembly_05.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa read_subsets/sample_06.fastq > assemblies/assembly_06.fasta

flye --nano-hq read_subsets/sample_07.fastq --threads "$threads" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07
./miniasm_and_minipolish.sh read_subsets/sample_08.fastq "$threads" > assemblies/assembly_08.gfa && any2fasta assemblies/assembly_08.gfa > assemblies/assembly_08.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa read_subsets/sample_09.fastq > assemblies/assembly_09.fasta

flye --nano-hq read_subsets/sample_10.fastq --threads "$threads" --out-dir assembly_10 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa assemblies/assembly_10.gfa && rm -r assembly_10
./miniasm_and_minipolish.sh read_subsets/sample_11.fastq "$threads" > assemblies/assembly_11.gfa && any2fasta assemblies/assembly_11.gfa > assemblies/assembly_11.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_12.gfa read_subsets/sample_12.fastq > assemblies/assembly_12.fasta

mv assemblies/assembly_2.fasta > assemblies/bad_assembly_2.fasta
mv assemblies/assembly_2.gfa > assemblies/bad_assembly_2.gfa
2 5 8 11 4 
2 4 '7' 11

trycycler cluster --assemblies assemblies/*.fasta --reads reads_qc/SRR8278838_sub.fastq.gz --out_dir trycycler --distance 0.02

mv trycycler/cluster_003 trycycler/bad_cluster_003

trycycler reconcile --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dir trycycler/cluster_001 --min_identity 97
mv ...
trycycler reconcile --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dir trycycler/cluster_002 
mv ...

trycycler msa --cluster_dir trycycler/cluster_001
trycycler msa --cluster_dir trycycler/cluster_002

trycycler partition --reads reads_qc/SRR8278838_sub.fastq.gz --cluster_dirs trycycler/cluster_*

bwa index trycycler.fasta
bwa mem -t 8 -a trycycler.fasta reads_qc/SRR8278185_1_sub.fastq.gz > alignments_1.sam
bwa mem -t 8 -a trycycler.fasta reads_qc/SRR8278185_2_sub.fastq.gz > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish trycycler.fasta filtered_1.sam filtered_2.sam > trycycler_polypolish.fasta
rm *.bwt *.pac *.ann *.amb *.sa *.sam


kraken2 --db kraken2/database/library/k2_standard_16gb_20230605 trycycler_polypolish.fasta  --report kraken.report --use-names


prokka --outdir prokka_results trycycler_polypolish.fasta 


kraken2/
reads/
reads_qc/
reads_subsets/
trycycles/

