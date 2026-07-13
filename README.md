# ntm-sgm-wgs

Scripts for benchmarked MAC genome recovery and lineage analysis.

Reads: NCBI BioProject `PRJNA1444780`, SRA study `SRP687912`.

Create the environment with `conda env create -f environment.yml`.

`scripts/run_recovery.sh` runs the two recovery routes and residual-mixture check. `scripts/near_mac_dilution.py` builds the dilution curves; `scripts/near_mac_two_route.py` runs the frozen 40-condition reconstruction and depth benchmark. The remaining scripts build the cohort, public context, lineage blocks and figures.

Software versions and accession tables are provided with the paper.
