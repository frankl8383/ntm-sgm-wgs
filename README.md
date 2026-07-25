# ntm-sgm-wgs

Code for MAC genome recovery, residual-mixture screening and lineage analysis.

Current release: `v2.6.2`.

Reads: NCBI BioProject `PRJNA1444780`; SRA study `SRP687912`.

```bash
conda env create -f environment.yml
conda activate ntm-sgm-wgs
make validate
```

Main stages:

1. `scripts/run_recovery.sh`: strict and meta recovery.
2. `scripts/residual_mixture.py`: residual within-MAC mixture screen.
3. `scripts/near_mac_dilution.py` and `scripts/near_mac_two_route.py`: mixture benchmarks.
4. `scripts/build_atlas.py`: public MAC atlas.
5. `scripts/accessory_families.py`, `scripts/define_blocks.py` and `scripts/validate_blocks.py`: accessory analysis.

Each script lists its inputs with `--help`. Frozen input manifests accompany the manuscript data archive.

`environment.yml` covers the Python code and tests. Install the external tools separately:

| Tool | Version |
|---|---:|
| fastp | 1.3.3 (local); 1.0.1 (external set) |
| Kraken2 / Bracken | 2.17.1 / 1.0.0 |
| Bowtie2 / samtools | 2.5.4 / 1.22.1 |
| SPAdes / QUAST | 4.2.0 / 5.3.0 |
| CheckM2 / GUNC | 1.1.0 / 1.0.6 |
| minimap2 / FastANI | 2.30 / 1.34 |
| NCBI Datasets | 18.33.1 |
| Prodigal / MMseqs2 / Bakta | 2.6.3 / 18.8cc5c / 1.11.2 |
| SKA2 / IQ-TREE | 0.5.1 / 3.0.1 |
| Mash / MashTree | 2.3 / 1.4.6 |

Reference databases:

| Use | Release | Source | Checksum |
|---|---|---|---|
| Kraken2 / Bracken | `k2_standard_16_GB_20260226` | [Standard-16 archive](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260226.tar.gz) | not retained |
| CheckM2 | downloaded 2026-07-10 | [Zenodo record 14897628](https://zenodo.org/records/14897628) | MD5 `07c10655620843b517d0df0c160d911f` |
| GUNC | ProGenomes 2.1, downloaded 2026-07-10 | [GUNC database](https://swifter.embl.de/~fullam/gunc/) | archive MD5 `bc93a855e0760aad5c4e5f2d0e26da46`; database MD5 `447c9330056b02f29f30fe81fe4af4eb` |
