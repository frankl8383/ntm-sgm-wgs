# ntm-sgm-wgs

Code for MAC genome recovery, residual-mixture screening and lineage analysis.

Reads: NCBI BioProject `PRJNA1444780`; SRA study `SRP687912`.

```bash
conda env create -f environment.yml
conda activate ntm-sgm-wgs
make validate
make figures
```

Analysis scripts are in `scripts/`. Figure code and source data are in `figures/`.
