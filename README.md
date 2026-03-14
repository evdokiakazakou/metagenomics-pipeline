# Metagenomics Pipeline

A pipeline for metagenomics data analysis using Oxford Nanopore and Illumina sequencing data.
Covers quality control, adapter trimming, taxonomic classification, genome assembly, binning,
and annotation.

---

## Table of Contents

1. [Download Public Data](#1-download-public-data)
2. [Quality Control](#2-quality-control)
3. [Adapter Trimming and Filtering](#3-adapter-trimming-and-filtering)
4. [Taxonomic Classification](#4-taxonomic-classification)
5. [Abundance Estimation with Bracken](#5-abundance-estimation-with-bracken)
6. [Genome Assembly with metaSPAdes](#6-genome-assembly-with-metaspades)
7. [Coverage Statistics](#7-coverage-statistics)
8. [Metagenome-Assembled Genomes with MetaBAT2](#8-metagenome-assembled-genomes-with-metabat2)
9. [Quality Assessment with CheckM](#9-quality-assessment-with-checkm)
10. [Phylogenetic Tree Visualisation](#10-phylogenetic-tree-visualisation)
11. [Annotation with Prokka](#11-annotation-with-prokka)

---

## Experimental Design

Before starting the analysis, document your experimental design in a table like the one below.
This helps keep track of which samples belong to which treatment group and which marker gene was used.

| Barcode   | Sample   | Treatment   | Marker | Platform |
|-----------|----------|-------------|--------|----------|
| barcode01 | sample_1 | treatment_A | ITS    | Nanopore |
| barcode02 | sample_2 | treatment_A | ITS    | Nanopore |
| barcode03 | sample_3 | treatment_B | 16S    | Nanopore |
| barcode04 | sample_4 | treatment_B | 16S    | Nanopore |

---

## 1. Download Public Data

Public metagenomics datasets are available through the NCBI Sequence Read Archive (SRA).

### 1.1 Find data with SRA Run Selector

Go to the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) and search
by accession number. For this course, use accession `PRJNA448333` from
[Cregger et al. 2018](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0618-5).

Note the individual sample accession numbers (e.g. SRR7301761) for the samples you want.

### 1.2 Download from the terminal

Use the SRA Toolkit to download data directly from the command line:

```bash
fastq-dump --split-files <accession_number>
```

- Replace `<accession_number>` with the actual sample accession (e.g. SRR7301761)
- `--split-files` splits paired-end reads into two separate files (_1 and _2)

---

## 2. Quality Control

FastQC checks the quality of raw sequencing reads and produces an HTML report per sample.
The report includes metrics such as per-base quality scores, GC content, and adapter content.

> Note: Nanopore data typically has lower quality scores than Illumina data.

### Nanopore — decompress and merge reads first

Each barcode folder may contain multiple .fastq.gz files. These need to be decompressed
and merged into a single file per barcode before running FastQC.

```bash
# Decompress all .fastq.gz files
# {} keeps the original filename, + processes each file one by one
find . -name '*.fastq.gz' -exec gunzip {} +

# Merge all reads per barcode into one file
for i in 01 02 03 04; do
    cat barcode${i}/*.fastq > barcode${i}_merged.fastq
done
```

### Run FastQC

```bash
module load gcc/14.2.0 fastqc/0.12.1

# For Nanopore merged files
fastqc *_merged.fastq

# For Illumina data
mkdir fastqc_out
for sample in *.fastq.gz; do
    fastqc ${sample} -o fastqc_out
done
```

Organise the output files:

```bash
mkdir fastqc
mv *.zip fastqc/
mv *.html fastqc/
```

---

## 3. Adapter Trimming and Filtering

Adapters are short sequences added during library preparation. They must be removed
before analysis because they are not part of the biological sequence.

### Nanopore — Porechop

Porechop finds and removes adapters from Nanopore reads. The `--discard_middle` flag
removes chimeric reads, which are reads that contain an adapter sequence in the middle,
indicating that two separate fragments were incorrectly joined.

```bash
module load gcc/14.2.0 miniconda3/24.7.1
source $CONDA_PROFILE/conda.sh

# Create and activate the conda environment (first time only)
conda create --name metagenomics_2026
conda activate metagenomics_2026

# Install Porechop (first time only)
conda install -c bioconda porechop

# Run Porechop on all barcodes
for i in 01 02 03 04; do
    porechop -i barcode${i}_merged.fastq \
             --discard_middle \
             -o barcode${i}_trimmed.fastq \
             --threads 4
done
```

Run FastQC again after trimming to confirm that quality has improved:

```bash
fastqc *_trimmed.fastq
```

### Nanopore — Custom filtering with NanoFilt

After Porechop, reads are filtered using NanoFilt, a tool designed specifically for
Oxford Nanopore data. It removes reads that are too short, too long, or have a low
average quality score. The thresholds are decided by observing the read length and
quality distribution plots from the previous step.

Recommended thresholds per marker gene:

| Marker | Min length | Max length | Min avg Q |
|--------|------------|------------|-----------|
| 16S    | 1300 bp    | 1600 bp    | 15        |
| ITS    | 200 bp     | 400 bp     | 15        |

```bash
# Install NanoFilt (first time only)
conda install -c bioconda nanofilt

mkdir filtered

# For 16S barcodes
for i in 03 04; do
    NanoFilt -l 1300 --maxlength 1600 -q 15 barcode${i}_trimmed.fastq \
    > ./filtered/barcode${i}_filtered.fastq
done

# For ITS barcodes
for i in 01 02; do
    NanoFilt -l 200 --maxlength 400 -q 15 barcode${i}_trimmed.fastq \
    > ./filtered/barcode${i}_filtered.fastq
done
```

| Flag | Description |
|------|-------------|
| `-l` | Minimum read length in bp |
| `--maxlength` | Maximum read length in bp |
| `-q` | Minimum average Phred quality score |

### Illumina — Trimmomatic

Trimmomatic trims adapter sequences and low-quality bases from Illumina paired-end reads.
It produces two output files per sample: paired (P, both reads survived) and unpaired
(U, only one read survived).

```bash
module load gcc/9.4.0-eewq4j6 trimmomatic/0.39-dlgljoz

trimmomatic PE -threads 4 -phred33 \
    sample_forward.fastq.gz sample_reverse.fastq.gz \
    sample_forward_P.fastq.gz sample_forward_U.fastq.gz \
    sample_reverse_P.fastq.gz sample_reverse_U.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:5 \
    TRAILING:5 \
    SLIDINGWINDOW:3:15 \
    MINLEN:100
```

| Parameter | Description |
|-----------|-------------|
| `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10` | Remove adapters using the TruSeq3 file; scan up to 30 nucleotides |
| `LEADING:5` | Trim 5 low-quality nucleotides from the start of each read |
| `TRAILING:5` | Trim 5 low-quality nucleotides from the end of each read |
| `SLIDINGWINDOW:3:15` | Scan 3 nucleotides at a time; discard if average quality drops below 15 |
| `MINLEN:100` | Discard reads shorter than 100 nucleotides after trimming |

Run FastQC again after trimming to confirm quality improvement.

---

## 4. Taxonomic Classification

Kraken2 classifies reads by comparing them against a reference database and assigning
each read to the most likely taxon.

### Install Kraken2

```bash
conda activate metagenomics_2026
conda install bioconda::kraken2
kraken2 --version
```

### Run Kraken2 for Nanopore 16S data (SILVA database)

```bash
for i in 03 04; do
    kraken2 --db ./kraken2_db/16S_database \
            --use-names \
            --confidence 0.1 \
            --report barcode${i}_kraken2_report.tabular \
            --output barcode${i}_kraken2_classification.tabular \
            ./filtered/barcode${i}_filtered.fastq
done
```

### Run Kraken2 for Nanopore ITS data (UNITE database)

```bash
for i in 01 02; do
    kraken2 --db ./kraken2_db/ITS_database \
            --use-names \
            --confidence 0.1 \
            --report barcode${i}_kraken2_report.tabular \
            --output barcode${i}_kraken2_classification.tabular \
            ./filtered/barcode${i}_filtered.fastq
done
```

### Run Kraken2 for Illumina data

```bash
kraken2 --db ./kraken2_db/standard_database \
        --paired \
        --use-names \
        --confidence 0.1 \
        --report sample_kraken2_report.tabular \
        --output sample_kraken2_classification.tabular \
        sample_forward_P.fastq.gz sample_reverse_P.fastq.gz
```

**Reference databases:**

| Marker   | Database | Version |
|----------|----------|---------|
| 16S rRNA | [SILVA](https://www.arb-silva.de/) | 138.2 |
| ITS      | [UNITE](https://unite.ut.ee/) | latest |
| Illumina | [Kraken2 Standard DB](https://benlangmead.github.io/aws-indexes/k2) | latest |

Visualise Kraken2 results interactively using [Krona](https://github.com/marbl/Krona).

---

## 5. Abundance Estimation with Bracken

Bracken uses Kraken2 classification reports to re-estimate the true abundance of each
taxon at a chosen taxonomic level.

### Install Bracken

```bash
conda install bioconda::bracken
```

### Run Bracken

```bash
for i in 01 02 03 04; do
    bracken -d ./kraken2_db/16S_database \
            -i barcode${i}_kraken2_report.tabular \
            -o barcode${i}_bracken.tabular \
            -r 1400 \
            -l G \
            -t 0
done
```

| Flag | Description |
|------|-------------|
| `-d` | Path to the Kraken2 database |
| `-i` | Input: Kraken2 report file |
| `-o` | Output: Bracken abundance file |
| `-r` | Read length (1400 bp for Nanopore) |
| `-l` | Taxonomic level: G=Genus, F=Family, S=Species, P=Phylum, C=Class, O=Order |
| `-t` | Minimum number of reads required for a taxon to be reported |

Run Bracken separately for each taxonomic level of interest and organise the results
into separate directories per level.

---

## 6. Genome Assembly with metaSPAdes

metaSPAdes assembles the trimmed paired-end reads into longer sequences called contigs.
A contig is a continuous DNA sequence reconstructed by overlapping reads.

```bash
metaspades.py \
    -1 sample_forward_P.fastq.gz \
    -2 sample_reverse_P.fastq.gz \
    -o metaspades_output
```

The main output is `metaspades_output/contigs.fasta`, which contains all assembled contigs.

---

## 7. Coverage Statistics

Before binning, we need to know how many reads map to each contig (coverage).
Higher coverage means more reads support that contig, which helps separate genomes.

```bash
module load gcc/14.2.0 bwa/0.7.17
module load gcc/14.2.0 samtools/1.19.2

# Index the assembled contigs
bwa index metaspades_output/contigs.fasta

# Map the original reads back to the contigs
bwa mem metaspades_output/contigs.fasta \
    sample_forward_P.fastq \
    sample_reverse_P.fastq > mapping.sam

# Convert SAM to BAM, sort, and index
samtools view -Sbu mapping.sam > mapping.bam
samtools sort mapping.bam -o mapping_sorted.bam
samtools index mapping_sorted.bam
```

---

## 8. Metagenome-Assembled Genomes with MetaBAT2

MetaBAT2 groups (bins) contigs into clusters that likely belong to the same organism,
producing metagenome-assembled genomes (MAGs). It uses both sequence composition
(tetranucleotide frequency) and coverage information.

```bash
metabat2 \
    -i metaspades_output/contigs.fasta \
    -o metabat2_output/bin \
    -m 1500 \
    --unbinned \
    --saveCls metabat2_output/bin.cls.tsv \
    --saveLog metabat2_output/metabat2.log
```

- `-m 1500` sets the minimum contig length to 1500 bp (shorter contigs are excluded)
- `--unbinned` saves contigs that could not be assigned to any bin

---

## 9. Quality Assessment with CheckM

CheckM estimates the completeness and contamination of each MAG by looking for
conserved marker genes that should be present in a complete genome.

```bash
checkm lineage_wf -t 4 -x fa metabat2_output/bin checkm_output
```

- `-t 4` uses 4 threads
- `-x fa` specifies the file extension of the bin files

---

## 10. Phylogenetic Tree Visualisation

CheckM produces a phylogenetic tree showing the evolutionary relationships between
the MAGs and reference genomes. The following tools can be used to visualise it:

| Tool | Type | Link |
|------|------|------|
| iTOL | Web-based (recommended) | http://itol.embl.de |
| FigTree | Desktop, newick format conversion | http://tree.bio.ed.ac.uk/software/figtree/ |
| Archaeopteryx | Desktop | https://sites.google.com/site/cmzmasek/home/software/archaeopteryx |
| EvolView | Web-based | https://www.evolgenius.info/evolview/ |

iTOL accepts newick format trees. Use FigTree to convert the tree format if needed
before uploading to iTOL.

---

## 11. Annotation with Prokka

Prokka annotates the MAGs by predicting genes and assigning functional information
to each predicted gene (e.g. gene name, product, database accession).

```bash
prokka \
    --outdir prokka_output \
    --prefix sample_name \
    metaspades_output/contigs.fasta
```

For broader functional annotation, [eggNOG-mapper](http://eggnog-mapper.embl.de/)
can be used as an alternative or in combination with Prokka.

---

## Acknowledgements

This pipeline is based on course material from the NGS and Metagenomics module of the Applied Bioinformatics MSc program (AUTH/IHU).
Original scripts developed by Dr. Spyros Papakostas (IHU, 2025-2026) and Dr. Loukas Theodosiou.
