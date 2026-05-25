# ntm-sgm-wgs

Code for the WGS reassessment workflow used in the associated SGM-NTM/MAC study.

## Contents

- `workflow/`: Snakemake workflow skeleton and modular rule files.
- `config/`: example paths, database settings and QC thresholds.
- `metadata/`: input templates and public-reference metadata templates.
- `scripts/`: analysis and figure scripts.

## Quick Start

Install Snakemake and mamba/conda, then prepare a samplesheet:

```bash
cp metadata/samplesheet.template.tsv metadata/samplesheet.tsv
```

Edit `metadata/samplesheet.tsv` and `config/config.yaml`, then run:

```bash
snakemake --use-conda --cores 8 --dry-run
snakemake --use-conda --cores 8 results/tables/preflight_samplesheet_check.tsv
```

## Data

Sequence reads are available from NCBI SRA. Accessions are listed in the paper and supplementary tables.
