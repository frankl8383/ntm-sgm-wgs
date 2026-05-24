# ntm-sgm-wgs

Minimal reproducibility code for a contamination-aware whole-genome reassessment workflow for presumed slowly growing nontuberculous mycobacterial clinical isolates.

This repository is intentionally code-only. It contains reusable workflow files, configuration examples, metadata templates and analysis scripts. It does not contain private data, generated outputs, large databases or author-only working files.

## What Is Included

- `workflow/`: Snakemake workflow skeleton and modular rule files.
- `workflow/envs/`: Conda environment files used by workflow modules.
- `config/`: Example configuration and QC threshold files.
- `scripts/`: Analysis and plotting scripts for sample validation, read-level taxonomy, assembly QC, species evidence integration, public-context ANI/phylogeny, AMR/stress review, and conservative mobile-element feature review.
- `metadata/`: Input templates and non-sensitive public-reference metadata templates.
- `CITATION.cff`: Citation metadata for this code repository.

## What Is Not Included

The repository deliberately excludes author-only working files, generated outputs, accession-management files, local paths, FASTQ files, genome assemblies, downloaded databases, public reference FASTA files and patient-level data.

## Data

Study data are managed through public sequence repositories and supplementary tables, not through this GitHub repository. Genome assemblies and large databases are regenerated or retrieved outside this code repository.

## Quick Start

Install Snakemake and mamba/conda, then create a project-specific samplesheet from the template:

```bash
cp metadata/samplesheet.template.tsv metadata/samplesheet.tsv
```

Edit `metadata/samplesheet.tsv` and `config/config.yaml` to point to local FASTQ files and local database paths. Then run a dry run:

```bash
snakemake --use-conda --cores 8 --dry-run
```

Run a preflight samplesheet check:

```bash
snakemake --use-conda --cores 8 results/tables/preflight_samplesheet_check.tsv
```

## Scope

This code supports the analysis framework used for contamination-aware WGS reassessment. It is not a clinical diagnostic workflow and should not be used to infer clinical antimicrobial susceptibility, transmission, complete plasmids, horizontal transfer, or formal taxonomic status without appropriate validation and supporting evidence.
