# ERSEM-ALLOC: reproducible standalone physiological experiments

**Deep S. Banerjee · Plymouth Marine Laboratory (PML)**
**Supported by the European Space Agency (ESA)**

This repository contains the Python workflows, configurations, reference outputs, observations and plotting materials used for the **standalone (zero-dimensional) PyFABM–ERSEM experiments** accompanying the manuscript *Cellular resource allocation reshapes phytoplankton physiology and seasonal ecosystem dynamics*.

ERSEM-ALLOC introduces a prognostic trade-off between effective carbon-fixation and nutrient-acquisition capacities while retaining the existing ERSEM nutrient-quota and chlorophyll-acclimation processes.

## Manuscript figures and experiments

- **Main Figure 3 — prescribed allocation:** balanced chlorophyll-to-carbon responses across irradiance, comparison with digitised *Thalassiosira pseudonana* observations from Geider et al. (1997), and P1 photosynthesis–irradiance responses for prescribed allocation fractions. The current figure-generation workflow includes a historically named output, `figures/supp_geider_pi.pdf`; this filename does **not** reflect the revised manuscript figure numbering.
- **Main Figure 4 — dynamic versus fixed allocation:** paired seven-day P1 experiments under nitrogen- and phosphorus-limited conditions, with constant incident shortwave radiation of 20 or 100 W m⁻² (`experiments/B4_dynamic_vs_fixed/`).
- **Reference outputs:** existing numerical results in `reference/`; outputs from subsequent runs in `results/` and experiment-specific result directories, as configured by the scripts.

## Repository contents

- `experiments/` — standalone experiment workflows, including the dynamic-versus-fixed comparison.
- `scripts/` — execution, plotting and release-manifest utilities.
- `configs/` — model/experiment configuration materials.
- `observations/` — observational inputs used by the standalone analyses.
- `reference/` — reference numerical outputs.
- `results/` — generated outputs and checks.
- `figures/` — figure products.
- `candidates/` — additional experimental analyses.
- `CITATION.cff` — citation metadata.

## Reproduce on Scylla

The experiments use PyFABM **2.1.4** with a compatible ERSEM-ALLOC installation and the `eat` Conda environment on Scylla.

From the repository root:

```bash
conda activate eat
bash scripts/run_all.sh
```

For the separate seven-day dynamic-versus-fixed experiment, choose a **new output directory** rather than reusing a reference-results directory:

```bash
B4_OUTPUT_DIR="$PWD/experiments/B4_dynamic_vs_fixed/results/B4_review_run" bash scripts/run_b4.sh
```

Use a separate output directory to preserve the supplied reference results.

## Licence

The original scripts are available under the BSD 3-Clause Licence (see `LICENSE`). ERSEM, FABM and observational datasets retain their respective terms.

## Credits and references

The experiments use the **European Regional Seas Ecosystem Model (ERSEM)** via **PyFABM**, the Python interface to the **Framework for Aquatic Biogeochemical Models (FABM)**. We acknowledge the ERSEM and FABM developer communities.

- **ERSEM:** Butenschön et al. (2016), *Geoscientific Model Development*. https://doi.org/10.5194/gmd-9-1293-2016
- **FABM:** Bruggeman and Bolding (2014), *Environmental Modelling & Software*. https://doi.org/10.1016/j.envsoft.2014.04.002
- **PyFABM documentation:** https://ersem.readthedocs.io/en/latest/tutorials/pyfabm.html
- **Geider observations:** Geider, MacIntyre and Kana (1997), *Marine Ecology Progress Series*. https://doi.org/10.3354/meps148187 (digitised Figure 3A observations).
- **ESA CAMP:** Climate Assessment and Model Parameterisation project. https://climate.esa.int/en/observations-and-modelling/camp/

See `CITATION.cff` for citation details.
