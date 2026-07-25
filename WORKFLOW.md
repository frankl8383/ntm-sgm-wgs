# Workflow

The main analysis stages are:

1. Build the public atlas with `scripts/build_atlas.py`.
2. Run strict and meta recovery with `scripts/run_recovery.sh`.
3. Summarize recovery gates with `scripts/summarize_recovery.py` and `scripts/build_cohort.py`; evaluate threshold multipliers with `scripts/threshold_sensitivity.py`.
4. Evaluate near-MAC mixtures with `scripts/near_mac_dilution.py`, `scripts/near_mac_two_route.py` and `scripts/random_pair_sensitivity.py`.
5. Summarize the frozen trees with `scripts/summarize_phylogeny.py`.
6. Build and test accessory families and intervals with `scripts/accessory_families.py`, `scripts/define_blocks.py` and `scripts/validate_blocks.py`.
7. Regenerate the figures with `make figures`.

Each script lists its required inputs with `--help`. Frozen source tables and figure inputs accompany the manuscript data archive.

`environment.yml` covers the Python scripts, tests and figures. Command-line tools were run in stage-specific environments at these versions:

| Tool | Version |
|---|---:|
| fastp | 1.3.3 (local); 1.0.1 (external evaluation) |
| Kraken2 / Bracken | 2.17.1 / 1.0.0 |
| Bowtie2 / samtools | 2.5.4 / 1.22.1 |
| SPAdes / QUAST | 4.2.0 / 5.3.0 |
| CheckM2 / GUNC | 1.1.0 / 1.0.6 |
| minimap2 / FastANI | 2.30 / 1.34 |
| NCBI Datasets | 18.33.1 |
| Prodigal / MMseqs2 / Bakta | 2.6.3 / 18.8cc5c / 1.11.2 |
| SKA2 / IQ-TREE | 0.5.1 / 3.0.1 |
| Mash / MashTree | 2.3 / 1.4.6 |

Database releases and checksums are listed in `DATABASES.md` and the supplementary version tables.
