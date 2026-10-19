# Original Section 2 B4 - seven-day dynamic vs fixed 0D allocation

Selected for **main manuscript Figure 3**, unlike the Geider/P-I supplementary physiology benchmark.

The original author's Section 2 PyFABM scripts (`code/`) and final manuscript Figure 3 plotting code (`figures/plot_fig3.py`) are included in this separate directory. The original source YAML is staged as `inputs/fabm.yaml` by `python scripts/bootstrap_scylla.py` in the parent PAPER repository. Nothing here modifies or submits the historic runs.

## Optional fresh numerical reproduction on Scylla

From `/work/dba/ERSEM-ALLOC-0D/PAPER` (with `eat` active), run `bash scripts/run_b4.sh`. This integrates the 7-day (dt=150 s) 20/100 W m-2, low-N/low-P/low-NP paired cases with dynamic/fixed fC. Output is **new** `experiments/B4_dynamic_vs_fixed/results/B4_dt150_fresh/` and plots to `figures/fig3.pdf` and `fig3.png`. This can take longer than the static Geider/P-I tests. It is intentionally separate from `scripts/run_all.sh`.

To preserve earlier reported CSV results for later comparison, run `python scripts/stage_b4_archived_results.py` (copy-only) from PAPER. This stages all `B4_*.csv` files under `reference/B4_dt150_original/`, not under fresh results. Compare the new results with the archived original prior to calling the B4 reproduction confirmed.

## Scientific scope

These are full-FABM Euler integrations with explicit fixed inorganic nutrient baths and host-updated P1 fC in dynamic arms. The fixed comparison holds P1 fC at 0.24 while leaving machinery-capacity factors active. It is not a separate Legacy compilation or a closed-budgets experiment. The full manuscript has additional 3D and satellite-observation analyses outside this 0D repository.
