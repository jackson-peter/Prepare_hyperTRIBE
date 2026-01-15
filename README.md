# Prepare_hyperTRIBE

This pipeline processes HyperTRIBE RNA-seq FASTQ (fq.gz) data into alignment, quantification, and per-base pileup outputs suitable for downstream analysis with the hyperTRIBER R package (for detecting differential RNA editing in HyperTRIBE experiments).​
STAR–Salmon–HyperTRIBE RNA-seq Pipeline

## Purpose

- Prepare HyperTRIBE RNA-seq data so it can be analyzed with hyperTRIBER.
- Perform:

    * Genome alignment with STAR, appropriate for edited reads.
    * Transcript quantification with Salmon.
    * Generation of mpileup-style base-count tables that can be reshaped into hyperTRIBER input (per-site coverage and base counts).​

## Inputs and assumptions

### Input data

- Paired-end RNA-seq FASTQ files (*.fq.gz) in the directory specified by fastq_dir.
Files share a common sample prefix and end with r1_suffix (e.g. _1P.fq.gz); the R2 suffix is inferred automatically.

- Reference files

- genome_fasta: Reference genome FASTA.

- genome_gtf: Corresponding annotation GTF for STAR.

- transcriptome_fasta: Transcriptome FASTA for Salmon indexing and quantification.​

## Workflow steps

- Indexing (if not already done)
    * Build STAR genome index from genome_fasta and genome_gtf.
    * Build Salmon index from transcriptome_fasta.
    * Index genome FASTA (samtools faidx) and generate .fai for downstream region splitting.

- Alignment and quantification

    * Align paired-end reads with STAR to produce coordinate-sorted BAM files plus STAR log files and splice junction tables.

    * Index BAM files with samtools.

    * Run salmon quant for each sample to obtain transcript-level expression estimates.​

    * HyperTRIBE-oriented mpileup

    * Extract chromosome names from one representative BAM.

    * Generate 1 Mbp genomic regions using the FASTA index.

    * Create a list of all BAM files.

    * For each region:

        + Run samtools mpileup on all BAMs.
        + Pipe to RNAeditR_mpileup2bases.pl to generate per-base counts text files, which can be converted into the tables used by hyperTRIBER for calling and testing editing sites.​

## Configuration

All key options are in config/config.yaml:

## Running the pipeline

```bash
# use snakemake 8.30
snakemake --use-conda --conda-frontend mamba -p --cores <N>
```

## Outputs:

    * STAR alignments and indices: results/alignment/star/.

    * Salmon quantification: results/quantification/salmon/.

    * Mpileup-derived per-region text files and a completion flag: results/variant/mpileup/, ready to be reshaped and imported into hyperTRIBER for differential editing analysis.