# Sub-Gaussian HDLSS PCA

Modular `amsart` manuscript and reproducibility package for arXiv submission.  This revision contains exactly 30 cited bibliography entries, expanding the literature coverage while leaving the theorem statements and numerical experiment outputs unchanged.

## Build

Requirements: TeX Live with `amsart`, `natbib`, `booktabs`, `mathtools`, `microtype`, `latexmk`, and BibTeX.

```bash
make
```

The compiled manuscript is `main.pdf`.

## Reproduce the numerical experiments

Python requirements are listed in `code/requirements.txt`.

```bash
make simulations
make
```

The checked-in results use master seed `20261109` and 1,000 replications per configuration.  `results/manifest.json` records the environment and design, while `results/simulation_summary.csv` contains means, Monte Carlo standard errors, medians, and 90th percentiles.

## Create the arXiv package

```bash
make arxiv
```

This creates `arxiv_submission.tar.gz` after a clean staging compile.  The paper uses a modular `sections/` layout and vector PDF figures.  The original manuscript appendices are intentionally not included.
