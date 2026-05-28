# Comparative Codon Usage Analysis Tools

This repository contains two pipelines for comparative codon usage analysis and codon optimization:

1. Chi-square codon bias pipeline
2. GLM-based codon bias pipeline

Both tools compare codon usage between target genes and selected reference datasets and support FASTA codon optimization.

---

# Features

- Comparative codon usage analysis
- Statistical comparison against reference datasets
- Scatterplot generation
- FASTA codon optimization
- Excel output generation

---

# Available Pipelines

## 1. Chi-square Pipeline

This pipeline uses contingency-table based chi-square statistics to compare synonymous codon distributions between target genes and reference datasets.

### Recommended for

- Biological interpretation
- Codon bias studies
- Comparative genomics
- Robust compositional analysis

### Statistics

- Chi-square test
- Benjamini-Hochberg multiple testing correction
- RSCU analysis
- log2(RSCU ratio)

---

## 2. GLM Pipeline

This pipeline uses generalized linear models (Poisson GLM) to model codon usage differences between target genes and references.

### Recommended for

- Model-based statistical analysis
- Effect size estimation
- Advanced statistical interpretation

### Statistics

- Poisson GLM
- Interaction modeling
- Multiple testing correction
- Codon effect coefficients

---

# Input Files

## Required

### Target Excel File

Excel file containing target genes/proteins.

Example:

| Targets |
|---|
| AOX1 |
| DAS1 |
| CAT1 |

The first column is used as input.

---

### FASTA Files

Required reference FASTA datasets:

- `op_cds_chr1-4.fasta`   #In this repo a mock dataset is used
- `NCBI_cds.fasta`

---

# References

Available reference datasets:

| Reference | Description |
|---|---|
| open | OpenPichia CDS dataset |
| ncbi | NCBI CDS dataset |
| both | Compare against both datasets |

---

# Using the GUI Applications

The pipelines are distributed as GUI executables.

To use the tools:

1. Double-click the `.exe` application
2. Select the target Excel file
3. Choose the reference dataset
4. Start the analysis

The program will automatically:

- perform codon usage analysis
- run statistical analyses
- generate plots
- export Excel result files
- create optimized FASTA files (if selected)

---

# Output

Each run creates a timestamped results directory:

```text
results_YYMMDD_HHMM/
```

Containing:

| File | Description |
|---|---|
| Excel results | Statistical output |
| Scatterplots | Codon/amino acid comparisons |
| Input | Copy of uploaded input files |
| Optimized FASTA | Codon-optimized sequences |

---
