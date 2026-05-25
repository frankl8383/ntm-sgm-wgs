# Workflow

Snakemake rules cover the core read QC, read-level taxonomy, assembly and genome-QC steps. Later ANI, species-evidence, AMR/stress and figure steps are provided as standalone scripts.

## Inputs

1. `metadata/samplesheet.tsv`
2. `config/config.yaml`
3. local databases listed in `config/databases.yaml`

## Basic commands

```bash
snakemake --use-conda --cores 8 --dry-run
snakemake --use-conda --cores 8 results/tables/preflight_samplesheet_check.tsv
snakemake --use-conda --cores 8 results/tables/read_qc_taxonomy_summary.tsv
snakemake --use-conda --cores 8 results/tables/assembly_qc_summary.tsv
snakemake --use-conda --cores 8 results/tables/final_wgs_species_reclassification_report_priority14.tsv
```

Several later public-context and figure-generation steps are script-driven because they were added after the core workflow skeleton. Representative examples include:

```bash
python scripts/build_integrated_ani_taxonomy_panel.py --help
python scripts/plot_figures_4_5_public_context.py --help
python scripts/review_ntm_resistance_loci.py --help
python scripts/summarize_genomad_results.py --help
```
