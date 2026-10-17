# Knee skeletal age acceleration and mortality in NHANES III

Code and derived data for:

> Sheng C, Zhang Y, Li X, Chen K. Deep learning–estimated skeletal age acceleration from knee radiographs and long-term mortality in older US adults: a cohort study. *Manuscript under review.*

We trained a deep learning model to estimate age from the bilateral knee radiographs of 4,875 NHANES III participants aged 60 years or older, and studied how knee skeletal age acceleration (KSAA), the part of predicted knee age not explained by chronological age, relates to mortality through 2019 and to 482 phenotypes.

![Study overview](docs/overview.png)

## Data

All inputs are public-use files from the National Center for Health Statistics:

- NHANES III (1988–1994) adult, examination and laboratory files: <https://wwwn.cdc.gov/nchs/nhanes/nhanes3/>
- NHANES III knee radiographs (4,908 films, about 171 GB) and the phase 2 radiograph reading file: <https://ftp.cdc.gov/pub/NHANES/XRays/Nhanes3/>
- Public-use Linked Mortality File, follow-up through 31 December 2019: <https://www.cdc.gov/nchs/data-linkage/mortality-public.htm>

The radiographs are not included here. Everything else needed for the statistical analyses is in `data/`:

| Folder | Content |
|---|---|
| `data/nhanes3/` | NHANES III files for the 4,908 participants with knee radiographs (DuckDB) |
| `data/predictions/` | Out-of-fold age predictions of the five backbones |
| `data/prepared/`, `data/images/`, `data/phewas/` | Variable extracts, image-processing logs and the PheWAS candidate pool |
| `data/derived/` | Analysis dataset, per-participant KSAA and PheWAS results, exactly as used in the manuscript |

See [`data/README.md`](data/README.md) and [`data/derived/README.md`](data/derived/README.md) for column definitions.

## Requirements

- R 4.5.2 (any R ≥ 4.3 should work) with:
  ```r
  install.packages(c("DBI", "duckdb", "ggplot2", "ggrepel", "here", "MASS", "mice", "mitools",
                     "patchwork", "png", "ragg", "survey", "survival"))
  ```
- Figure 4A only: Python ≤ 3.10 with CLARITE 2.3.6, PyMuPDF and Pillow:
  ```bash
  conda create -n clarite python=3.9 -y
  conda run -n clarite pip install clarite==2.3.6 PyMuPDF Pillow
  ```
- Image processing and model training only: Python 3 with PyTorch 2.1 and an NVIDIA GPU (we used one RTX 3090).

Exact versions are listed in [`environment/`](environment/). The analyses were run on macOS.

## Reproducing the results

From the repository root:

```bash
CLARITE_PY=$(conda run -n clarite which python) bash run_all.sh
```

This runs the scripts in `code/analysis/` in numbered order (01–20) and takes about 5 minutes on a laptop. Tables are written to `results/analysis/`, figures to `results/figures/` and console logs to `results/logs/`. The last step, `verify_derived.R`, compares the regenerated results with `data/derived/` and stops with an error if anything differs; a successful run ends with `ALL IDENTICAL`.

To re-run selected scripts after a full run, give their numbers, e.g. `bash run_all.sh 13 14`.

| Manuscript item | Output (in `results/`) | Script |
|---|---|---|
| Table 1, Tables S1 and S3 | `analysis/tables/table1.html`, `tableS1.html`, `tableS3.html` | `13_manuscript_tables.R` |
| Table 2 | `analysis/weighted/table2.csv` | `06_weighted.R` |
| Tables S4, S5, S7, S10 | `analysis/weighted/tableS*.csv` | `06_weighted.R` |
| Table S6 | `analysis/weighted/tableS6.csv` | `07_weighted_boot.R` |
| Tables S8 and S9 | `analysis/tables/tableS8.html`, `tableS9.html` | `08`–`11_phewas*.R` |
| Figure 1 | `figures/Figure1/Figure1.pdf` | `15_fig1.R` |
| Figure 2, Fig. S3 | `figures/Figure2.pdf`, `FigureS3.pdf` | `16_fig2_mortality.R` |
| Figure 3 | `figures/Figure3.pdf` | `17_fig3_incremental.R` |
| Figure 4 | `figures/Figure4.pdf` | `18_fig4_phewas.R`, `19_fig4_phewas_clarite.py` |
| Fig. S1 | `figures/FigureS1.pdf` | `20_figS1_biascorr.R` |
| Fig. S2 | `figures/FigureS2.pdf` (needs the radiographs and a GPU) | `code/interpretability/` |
| Other numbers in the text | `logs/*.txt` | `02`, `12`, `14` |

## Knee age model

The model adapts [Mammo-AGE](https://doi.org/10.1038/s41467-025-65923-5) to read both knees together. The steps below produced the files in `data/predictions/` and are not needed to reproduce the analyses.

1. `code/images/`: download the radiographs, detect both knees with YOLO11 and crop regions of interest ([README](code/images/README.md)).
2. `code/model/train_5fold.sh`: train the five backbones with participant-level five-fold cross-fitting, about 100 GPU-hours ([README](code/model/README.md)).
3. `code/analysis/02_ksaa_ensemble.R`: combine the backbones and correct for age bias to obtain KSAA.

Models trained on the full cohort (`full_<backbone>.pt`, 559 MB) will be deposited on Zenodo (DOI to be added); [`models/README.md`](models/README.md) describes how to apply them to new radiographs. Performance on weight-bearing or digital radiographs and in adults younger than 60 years has not been evaluated.

## Citation

If you use this code or data, please cite the article above (citation will be updated on publication).

## License

The code is released under the [MIT License](LICENSE). Files in `code/model/` adapted from Mammo-AGE remain under the Apache License 2.0 ([`LICENSE_Mammo-AGE`](code/model/LICENSE_Mammo-AGE)). NHANES data are public-use data from NCHS and are subject to the NCHS data use restrictions.

## Contact

Xiangchun Li (lixiangchun@tmu.edu.cn) or Kexin Chen (chenkexin@tmu.edu.cn), Tianjin Medical University Cancer Institute and Hospital. Questions about the code can also be posted as GitHub issues.
