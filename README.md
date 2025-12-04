# motif_rater

A utility for summarizing motif content in genome FASTA files.

## Usage

```bash
cargo run -- \
  --genome-fasta-file genome1.fa --genome-fasta-file genome2.fa \
  --motif ATGC --motif-file motifs.txt
```

Inputs mirror common CoverM FASTA options: provide genome files directly, supply a file with newline-delimited paths via `--genome-fasta-list`, or scan directories with `--genome-fasta-directory`.

Output columns (tab-separated) per genome/motif combination:

1. `genome`: FASTA file name.
2. `length`: total number of bases.
3. `gc_content`: GC fraction across the genome.
4. `motif`: motif label.
5. `count`: occurrences of the motif plus its reverse complement.
6. `expected`: expected occurrences given GC content and genome length.
7. `observed_rate`: count divided by possible windows.
8. `poisson_p_value`: upper-tail Poisson probability of observing at least `count` given `expected`.
