# Reproducing the analysis workflow

This project contains a Snakemake-based analysis and supporting scripts for contamination-aware WGS reassessment of presumed SGM-NTM isolates.

## Inputs

1. `metadata/samplesheet.tsv`: local FASTQ manifest created from `metadata/samplesheet.template.tsv`.
2. `config/config.yaml` and `config/qc_thresholds.yaml`: paths and thresholds.
3. Local database installations matching `config/databases.yaml`.
4. Public genome metadata and FASTA files downloaded locally by the user.

## Minimal command sequence

The project was built incrementally. For a fresh reproduction on a machine with the databases available, start with:

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

## Interpretation Constraints

Public-context phylogenies should not be interpreted as transmission analyses. Genome-only AMR/stress feature review is not clinical resistance prediction. Conservative mobile-element-associated feature review from short-read draft assemblies does not prove complete plasmids, horizontal transfer or epidemiological linkage.
