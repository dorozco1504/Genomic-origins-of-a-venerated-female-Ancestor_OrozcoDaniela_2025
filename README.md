# Ancient DNA Preprocessing and Variant Calling Pipeline

This repository contains a **Snakemake workflow** used for processing ancient DNA (aDNA) data. It was developed as part of the analyses presented in the study:

> **Genomic insights on the worship of a female Ancestor in ancient Mexico**  
> *Daniela Orozco Pérez et al., 2025*

The pipeline includes steps for adapter trimming, alignment to a reference genome, quality filtering, variant calling with ANGSD, and conversion to PLINK format for downstream population genetics analyses.

---

## 📦 Features

- Adapter trimming with `AdapterRemoval`
- Alignment with `BWA aln` to a human reference genome
- Quality filtering and duplicate removal using `samtools`
- SNP calling at predefined sites using `ANGSD`
- Conversion of results to PLINK binary format for downstream analyses

---

## 🧪 Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/) ≥ 7.0  
- [Conda](https://docs.conda.io/en/latest/) for environment management

> All dependencies are provided through the included Conda environment file (`environment.yml`).

---

## 📁 Repository Structure

```
.
├── Snakefile               # Main workflow
├── config.yaml             # User configuration (paths, threads, dataset info)
├── environment.yml         # Conda environment definition
└── results/                # Output files (set via config.yaml)
```
---

## ⚙️ Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/dorozco1504/aDNA-pipeline.git
   cd aDNA-pipeline
   ```

2. **Create the environment**:
   ```bash
   conda env create -f environment.yml --name adna_env
   conda activate adna_env
   ```

3. **Edit `config.yaml`** to reflect your input/output directories and dataset parameters.

4. **Run the workflow**:
   ```bash
   snakemake --cores 16 --use-conda
   ```

---

## 🔧 Configuration

Minimal example for `config.yaml`:

```yaml
rounds:
  - CdV

threads:
  adapter_removal: 8
  bwa: 8
  samtools: 8

Paths:
  data_dir: /path/to/raw_fastq
  output_dir: /path/to/output
  sites: /path/to/predefined_sites.bed
```

---

## 🧬 Output Files

- Collapsed FASTQs from AdapterRemoval
- BAMs at different filtering stages
- Individual-level alignment statistics
- Haplotype calls from ANGSD
- `.bed`, `.bim`, and `.fam` files for use with PLINK

---

## 📚 Citing This Pipeline

If you use this workflow, please cite the original study:

**Orozco Pérez, D.**, et al. (2025). *Genomic insights on the worship of a female Ancestor in ancient Mexico*.  
[Journal Name, Volume(Issue), Pages]. [DOI or link if available]

Also cite the software used: AdapterRemoval, BWA, Samtools, ANGSD, and PLINK.

---
