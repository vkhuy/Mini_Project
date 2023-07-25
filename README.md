# MINI PROJECT - MICROBIAL GENOME ANALYSIS

## Introduction
This repository is a bioinformatics analysis for microbial genome assembly, annotation and taxonomic classification. The pipeline is able to assemble short reads and long reads (hybrid assembly).

## Installation and setup

Because there are many dependencies to install across all of the different analytical tools the first time you run the tool it can take some time to install everything via conda. Because of this I have added a script to install environments for the six modules `setup_environment.sh`.

```
# Clone the git repository and enter it
git clone https://github.com/vkhuy/Mini_Project.git
cd Mini_project

# Setup and installation - this may take some time
bash setup_environments.sh
```

## Usage

Run the script by providing the two input FASTQ files as arguments:
```
bash script.sh input_file1.fastq.gz input_file2.fastq.gz input_file3.fastq.gz 
```
- `input_file1.fastq.gz` and `input_file2.fastq.gz` are Illumina paired end reads.
- `input_file3.fastq.gz` is ONT reads.

## Modules/Workflows

This pipeline has six main modules:

1. Quality control
2. Assembly
3. Long read polishing
4. Short read polishing
5. Taxonomic classification
6. Annotation

### 1. Quality control

The goal of read QC is to discard low-quality reads and/or trim off low-quality regions of reads. This will make them easier to use in later steps (assembly and polishing).

For Illumina read QC, use [fastp](https://github.com/OpenGene/fastp) to remove adapters and trim off low-quality bases.

For ONT read QC, use [Filtlong](https://github.com/rrwick/Filtlong) to remove short reads to improve N50 and throw out the worst 10% of reads. After these QC steps, you should be left with an ONT read set with a much better N50 but still plenty of depth. Then it will visiualize ONT read QC with [FastQC](https://github.com/s-andrews/FastQC).

### 2. Assembly 

In this pipeline, I recommend use [Trycycler](https://github.com/rrwick/Trycycler). Trycycler is a tool that takes as input multiple separate long-read assemblies of the same genome (e.g. from different assemblers or different read subsets) and produces a consensus long-read assembly.

In brief, Trycycler does the following:

- Clusters the contig sequences, so the user can distinguish complete contigs (i.e. those that correspond to an entire replicon) from spurious and/or incomplete contigs.
- Reconciles the alternative contig sequences with each other and repairs circularisation issues.
- Performs a multiple sequence alignment (MSA) of the alternative sequences.
- Constructs a consensus sequence from the MSA by choosing between variants where the sequences differ.

There are 6 steps to assembly:
- [Generating assemblies](https://github.com/rrwick/Trycycler/wiki/Generating-assemblies): In this step, there are 3 assemblers: [Flye](https://github.com/fenderglass/Flye), [Miniasm+Minipolish](https://github.com/rrwick/Minipolish), [Raven](https://github.com/lbcb-sci/raven).
- [Clustering contigs](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs)
- [Reconciling contigs](https://github.com/rrwick/Trycycler/wiki/Reconciling-contigs)
- [Multiple sequence alignment](https://github.com/rrwick/Trycycler/wiki/Multiple-sequence-alignment)
- [Partitioning reads](https://github.com/rrwick/Trycycler/wiki/Partitioning-reads)
- [Generating a consensus](https://github.com/rrwick/Trycycler/wiki/Generating-a-consensus)

For these commands to complete, it may be necessary to delete or repair some of the assembly/cluster sequences by manually inspect. If you have a different sample of my sample, you should read the step 1, 2, 3 above to fix it.

Sometime, [Flye assembly produces differentiating output on the same input](https://github.com/fenderglass/Flye/issues/509). It is normal, just do assembly with Flye again, you will have the same output like me.

```
flye --nano-hq read_subsets/sample_"$sample".fastq --threads "$threads" --out-dir assembly_"$sample" && cp assembly_"$sample"/assembly.fasta assemblies/assembly_"$sample".fasta && cp assembly_"$sample"/assembly_graph.gfa assemblies/assembly_"$sample".gfa && rm -r assembly_"$sample"
```

### 3. Long read polishing

[Medaka](https://github.com/nanoporetech/medaka) polishing aims to fix as many remaining errors in your assembly using only ONT reads. The goal is to produce the best possible ONT-only assembly.

### 4. Short read polishing

The goal of short-read polishing is to fix these small-scale errors using Illumina reads. [Polypolish](https://github.com/rrwick/Polypolish) is a very safe polisher to do it.

First, align Illumina reads (separately for the `_1` and `_2` files) to your assembly using the all-alignments-per-read option. Then, Polypolish will fix the errors on the draft assembly.

### 5. Taxonomic classification

After obtaining the complete assembly sequence, the target species will be identified using the [kraken2](https://github.com/DerrickWood/kraken2) with a 16GB database including reference genomes of bacteria, viruses, humans, etc of NCBI. This database will dowload automatic in the scripts so be patient and wait until it done.

### 6. Annotation

The resulting bacterial assembly is furthermore annotated using [Prokka](https://github.com/tseemann/prokka).

## Output directory structure and file naming

For more details, refer to the file table section in an HTML report generated by the pipeline. Make sure you have the similar output.

```
Mini_Project
├── Modules
│   ├── Module01_quality_control.sh
│   ├── Module02_assembly.sh
│   ├── Module03_long_read_polishing.sh
│   ├── Module04_short_read_polishing.sh
│   ├── Module05_taxonomic_classification.sh
│   └── Module06_Annotation.sh
├── any2fasta
├── assembled
│   ├── trycycler.fasta
│   ├── trycycler_medaka.fasta
│   └── trycycler_medaka_polypolish.fasta
├── assemblies
│   ├── assembly_01.fasta
│   ├── assembly_01.gfa
│   ├── assembly_03.fasta
│   ├── assembly_03.gfa
│   ├── assembly_05.fasta
│   ├── assembly_05.gfa
│   ├── assembly_06.fasta
│   ├── assembly_06.gfa
│   ├── assembly_08.fasta
│   ├── assembly_08.gfa
│   ├── assembly_09.fasta
│   ├── assembly_09.gfa
│   ├── assembly_10.fasta
│   ├── assembly_10.gfa
│   ├── assembly_12.fasta
│   ├── assembly_12.gfa
│   └── bad_assembly
│       ├── assembly_02.fasta
│       ├── assembly_02.gfa
│       ├── assembly_04.fasta
│       ├── assembly_04.gfa
│       ├── assembly_07.fasta
│       ├── assembly_07.gfa
│       ├── assembly_11.fasta
│       └── assembly_11.gfa
├── figtree.jar
├── kraken2
│   └── k2_standard_16gb_20230605
│       ├── database100mers.kmer_distrib
│       ├── database150mers.kmer_distrib
│       ├── database200mers.kmer_distrib
│       ├── database250mers.kmer_distrib
│       ├── database300mers.kmer_distrib
│       ├── database50mers.kmer_distrib
│       ├── database75mers.kmer_distrib
│       ├── hash.k2d
│       ├── inspect.txt
│       ├── ktaxonomy.tsv
│       ├── opts.k2d
│       ├── seqid2taxid.map
│       └── taxo.k2d
├── kraken2_report
│   ├── SRR8278838_sub.kraken
│   └── kraken.report
├── miniasm_and_minipolish.sh
├── prokka_results
│   ├── PROKKA_07242023.err
│   ├── PROKKA_07242023.faa
│   ├── PROKKA_07242023.ffn
│   ├── PROKKA_07242023.fna
│   ├── PROKKA_07242023.fsa
│   ├── PROKKA_07242023.gbf-r
│   ├── PROKKA_07242023.gbk
│   ├── PROKKA_07242023.gff
│   ├── PROKKA_07242023.log
│   ├── PROKKA_07242023.sqn
│   ├── PROKKA_07242023.tbl
│   ├── PROKKA_07242023.tsv
│   └── PROKKA_07242023.txt
├── qc_report
│   ├── SRR8278838_sub_fastqc.html
│   ├── SRR8278838_sub_fastqc.zip
│   ├── fastp.html
│   └── fastp.json
├── read_subsets
│   ├── sample_01.fastq
│   ├── sample_02.fastq
│   ├── sample_03.fastq
│   ├── sample_04.fastq
│   ├── sample_05.fastq
│   ├── sample_06.fastq
│   ├── sample_07.fastq
│   ├── sample_08.fastq
│   ├── sample_09.fastq
│   ├── sample_10.fastq
│   ├── sample_11.fastq
│   └── sample_12.fastq
├── reads
│   ├── SRR8278185_1_sub.fastq.gz
│   ├── SRR8278185_2_sub.fastq.gz
│   └── SRR8278838_sub.fastq.gz
├── reads_qc
│   ├── SRR8278185_1_sub.fastq.gz
│   ├── SRR8278185_2_sub.fastq.gz
│   ├── SRR8278838_sub.fastq.gz
│   └── unpair.fastq.gz
├── script.sh
├── setup_environment.sh
└── trycycler
    ├── bad_cluster_003
    │   └── 1_contigs
    │       ├── A_contig_2.fasta
    │       └── G_contig_1.fasta
    ├── cluster_001
    │   ├── 1_contigs
    │   │   ├── A_contig_1.fasta
    │   │   ├── B_Utg1910.fasta.bad
    │   │   ├── C_utg000001c.fasta
    │   │   ├── D_Utg1850.fasta
    │   │   ├── E_utg000001c.fasta
    │   │   ├── F_Utg1744.fasta
    │   │   ├── G_contig_3.fasta
    │   │   └── H_Utg1908.fasta
    │   ├── 2_all_seqs.fasta
    │   ├── 3_msa.fasta
    │   ├── 4_reads.fastq
    │   ├── 5_chunked_sequence.gfa
    │   ├── 6_initial_consensus.fasta
    │   ├── 7_final_consensus.fasta
    │   └── 8_medaka.fasta
    ├── cluster_002
    │   ├── 1_contigs
    │   │   ├── A_contig_3.fasta.bad
    │   │   ├── B_Utg1908.fasta
    │   │   ├── C_utg000002c.fasta
    │   │   ├── D_Utg1848.fasta
    │   │   ├── E_utg000002c.fasta
    │   │   ├── F_Utg1742.fasta.bad
    │   │   ├── G_contig_2.fasta.bad
    │   │   └── H_Utg1906.fasta
    │   ├── 2_all_seqs.fasta
    │   ├── 3_msa.fasta
    │   ├── 4_reads.fastq
    │   ├── 5_chunked_sequence.gfa
    │   ├── 6_initial_consensus.fasta
    │   ├── 7_final_consensus.fasta
    │   └── 8_medaka.fasta
    ├── clusters.png
    ├── contigs.newick
    └── contigs.phylip
```
