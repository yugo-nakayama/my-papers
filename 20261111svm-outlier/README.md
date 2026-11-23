# Outliers and the Hard-Margin Linear SVM — arXiv source

[![arXiv](https://img.shields.io/badge/arXiv-2611.13869-b31b1b.svg)](https://arxiv.org/abs/2611.13869)

- **arXiv**: [arXiv:2611.13869](https://arxiv.org/abs/2611.13869) [stat.ML]
- **Title**: Outliers and the Hard-Margin Linear SVM in High-Dimension, Low-Sample-Size Settings
- **Author**: Yugo Nakayama

---

This directory is a self-contained manuscript project following the modular layout of `yugo-nakayama/my-papers`.

## Build

```bash
make
```

Requirements: TeX Live with `amsart`, `amsmath`, `amssymb`, `amsthm`, `mathtools`, `natbib`, `geometry`, `booktabs`, `graphicx`, `caption`, `subcaption`, `enumitem`, `microtype`, and `hyperref`.

## Reproduce numerical experiments

```bash
python3 -m pip install -r experiments/requirements.txt
make sim          # manuscript dimensions and replication counts
make sim-quick    # smoke test with smaller dimensions and fewer repetitions
```

The fixed default seed is `20261029`. The script writes raw CSV files, generated LaTeX tables, figures, and `results/run_metadata.json`.

## Package for arXiv

```bash
make arxiv
```

This creates `arxiv_submission.tar.gz` containing `main.tex`, the modular section files, bibliography, PDF figures, generated table fragments, experiment code, raw numerical outputs, and run metadata. Build products and the source PDF used for reconstruction are intentionally excluded.

## File map

- `main.tex`: document class, packages, theorem environments, macros, title, and section inputs.
- `sections/`: manuscript sections.
- `figures/`: generated PDF figures with caption-safe whitespace.
- `experiments/run_experiments.py`: all seven studies.
- `results/`: CSV values, LaTeX table fragments, software/seed metadata.
- `references.bib`: BibTeX database using author-year citations.
