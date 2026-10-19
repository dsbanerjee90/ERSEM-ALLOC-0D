(base) deepsankarbanerjee@Deeps-MacBook-Pro Downloads $ cat README.md
# ERSEM-ALLOC: reproducible standalone physiological experiments

**Deep S. Banerjee · Plymouth Marine Laboratory (PML)**
**Supported by the European Space Agency (ESA)**

This repository contains the Python workflows, configurations, reference outputs, observations and plotting materials used for the **standalone (zero-dimensional) PyFABM–ERSEM experiments** accompanying the manuscript *Cellular resource allocation reshapes phytoplankton physiology and seasonal ecosystem dynamics*.

ERSEM-ALLOC introduces a prognostic trade-off between effective carbon-fixation and nutrient-acquisition capacities while retaining the existing ERSEM nutrient-quota and chlorophyll-acclimation processes.

## Manuscript figures and experiments

- **Main Figure 3 — prescribed allocation:** balanced chlorophyll-to-carbon responses across irradiance, comparison with digitised *Thalassiosira pseudonana* observations from Geider et al. (1997), and P1 photosynthesis–irradiance responses for prescribed allocation fractions. The current figure-generation workflow includes a historically named output, `figures/supp_geider_pi.pdf`; this filename does **not** reflect the revised manuscript figure numbering.
- **Main Figure 4 — dynamic versus fixed allocation:** paired seven-day P1 experiments under nitrogen- and phosphorus-limited conditions, with constant incident shortwave radiation of 20 or 100 W m⁻² (`experiments/B4_dynamic_vs_fixed/`).
- **Reference outputs:** existing numerical results in `reference/`; outputs from subsequent runs in `results/` and experiment-specific result directories, as configured by the scripts.

The experiment outputs and final manuscript figure files should be checked against the archived figure-generation scripts before a public release.

## Repository contents

- `experiments/` — standalone experiment workflows, including the dynamic-versus-fixed comparison.
- `scripts/` — execution, plotting and release-manifest utilities.
- `configs/` — model/experiment configuration materials.
- `observations/` — observational inputs used by the standalone analyses.
- `reference/` — reference numerical outputs.
- `results/` — generated outputs and checks.
- `figures/` — figure products.
- `provenance/` — available provenance records.
- `candidates/` — candidate and intermediate materials; inclusion in a release requires review.
- `CITATION.cff` — citation metadata (to be checked before tagging a release).
- `LICENSE_DECISION_REQUIRED.md` — pending licensing decision.

## Reproduce on Scylla

These workflows currently depend on an **existing ERSEM-ALLOC-compatible ERSEM + FABM/PyFABM installation** and the local `eat` Conda environment. The manuscript reports PyFABM version **2.1.4**; the underlying ERSEM/FABM source revision and executable-library provenance must be recorded before describing the environment as independently reproducible.

From the repository root:

```bash
conda activate eat
bash scripts/run_all.sh
python scripts/release_manifest.py
```

For the separate seven-day dynamic-versus-fixed experiment, choose a **new output directory** rather than reusing a reference-results directory:

```bash
B4_OUTPUT_DIR="$PWD/experiments/B4_dynamic_vs_fixed/results/B4_review_run" bash scripts/run_b4.sh
```

Inspect each script's output paths before execution. Preserve the reference outputs and avoid overwriting existing runs. The commands above document the current Scylla workflow; an installation-independent test has **not** yet been established.

## Scope and source-code availability

**This directory is not yet a complete release of the ERSEM-ALLOC model.** In particular, it does not contain the modified ERSEM Fortran implementation, the NEMO–FABM coupling changes, or the full three-dimensional regional simulation and evaluation workflow. Consequently, the standalone experiments cannot be rerun from this repository alone without access to the matching modified model installation.

Before manuscript submission, the authors intend to identify and archive the exact modified source files or an appropriate patch against a documented upstream revision, together with build instructions, model configuration and version provenance. Access for editors and reviewers will be arranged through an authorised private repository or secure review package. Public release and DOI registration will follow resolution of licensing, redistribution permissions and repository validation.

## Licensing and release status

**Pre-release; licensing pending.** See `LICENSE_DECISION_REQUIRED.md`. Do not assume that the original ERSEM/FABM code, modified source files or third-party observational datasets can be redistributed under a newly chosen licence without checking their respective terms. Exclude credentials, institution-only data and machine-specific private paths from a public archive.

A private GitHub repository can be used for preparation and controlled review. **Zenodo's automated GitHub archiving works with public repositories, not private ones.** After the repository is cleared for public release, connect it to Zenodo and create a versioned GitHub release. Add the resulting archived DOI to `CITATION.cff`, the README and the manuscript's code-availability statement.

## Credits and references

The experiments use the **European Regional Seas Ecosystem Model (ERSEM)** via **PyFABM**, the Python interface to the **Framework for Aquatic Biogeochemical Models (FABM)**. We acknowledge the ERSEM and FABM developer communities.

- **ERSEM:** Butenschön et al. (2016), *Geoscientific Model Development*. https://doi.org/10.5194/gmd-9-1293-2016
- **FABM:** Bruggeman and Bolding (2014), *Environmental Modelling & Software*. https://doi.org/10.1016/j.envsoft.2014.04.002
- **PyFABM documentation:** https://ersem.readthedocs.io/en/latest/tutorials/pyfabm.html
- **Geider observations:** Geider, MacIntyre and Kana (1997), *Marine Ecology Progress Series*. https://doi.org/10.3354/meps148187 (digitised Figure 3A observations).
- **ESA CAMP:** Climate Assessment and Model Parameterisation project. https://climate.esa.int/en/observations-and-modelling/camp/

Once an authorised versioned archive is available, cite the release DOI alongside the original ERSEM/FABM and observational references. See `CITATION.cff` for repository citation metadata.
