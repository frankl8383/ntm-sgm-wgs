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

Each script lists its inputs with `--help`. Software and database versions are reported in the supplementary tables.
