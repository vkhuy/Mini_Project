#!usr/bin/env bash

# Annotation with prokka
mamba run -n prokka prokka --outdir prokka_results --genus helicobacter --species pylori  trycycler_medaka_polypolish.fasta
