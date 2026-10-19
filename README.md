# Knee skeletal age acceleration and mortality in NHANES III

Code and derived data for:

> Sheng C, Zhang Y, Li X, Chen K. Deep learning–estimated skeletal age acceleration from knee radiographs and long-term mortality in older US adults: a cohort study. *Manuscript under review.*

![Study overview](docs/overview.png)

## Data

All inputs are public-use files from the National Center for Health Statistics:

- NHANES III (1988–1994) adult, examination and laboratory files: <https://wwwn.cdc.gov/nchs/nhanes/nhanes3/>
- NHANES III knee radiographs (4,908 films, about 171 GB) and the phase 2 radiograph reading file: <https://ftp.cdc.gov/pub/NHANES/XRays/Nhanes3/>
- Public-use Linked Mortality File, follow-up through 31 December 2019: <https://www.cdc.gov/nchs/data-linkage/mortality-public.htm>

The data and the analysis files derived from them (`data/`) are not included in this repository.

## Requirements

- R 4.5.2 with:
  
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

With the input files in `data/`, run from the repository root:

```bash
CLARITE_PY=$(conda run -n clarite which python) bash run_all.sh
```

This runs the scripts in `code/analysis/` in numbered order (01–20) and takes about 5 minutes on a laptop. Tables are written to `results/analysis/`, figures to `results/figures/` and console logs to `results/logs/`. The last step, `verify_derived.R`, compares the regenerated results with `data/derived/` and stops with an error if anything differs; a successful run ends with `ALL IDENTICAL`.

To re-run selected scripts after a full run, give their numbers, e.g. `bash run_all.sh 13 14`.

| Manuscript item           | Output (in `results/`)                                        | Script                                          |
| ------------------------- | ------------------------------------------------------------- | ----------------------------------------------- |
| Table 1, Tables S1 and S3 | `analysis/tables/table1.html`, `tableS1.html`, `tableS3.html` | `13_manuscript_tables.R`                        |
| Table 2                   | `analysis/weighted/table2.csv`                                | `06_weighted.R`                                 |
| Tables S4, S5, S7, S10    | `analysis/weighted/tableS*.csv`                               | `06_weighted.R`                                 |
| Table S6                  | `analysis/weighted/tableS6.csv`                               | `07_weighted_boot.R`                            |
| Tables S8 and S9          | `analysis/tables/tableS8.html`, `tableS9.html`                | `08`–`11_phewas*.R`                             |
| Figure 1                  | `figures/Figure1/Figure1.pdf`                                 | `15_fig1.R`                                     |
| Figure 2, Fig. S3         | `figures/Figure2.pdf`, `FigureS3.pdf`                         | `16_fig2_mortality.R`                           |
| Figure 3                  | `figures/Figure3.pdf`                                         | `17_fig3_incremental.R`                         |
| Figure 4                  | `figures/Figure4.pdf`                                         | `18_fig4_phewas.R`, `19_fig4_phewas_clarite.py` |
| Fig. S1                   | `figures/FigureS1.pdf`                                        | `20_figS1_biascorr.R`                           |
| Fig. S2                   | `figures/FigureS2.pdf` (needs the radiographs and a GPU)      | `code/interpretability/`                        |
| Other numbers in the text | `logs/*.txt`                                                  | `02`, `12`, `14`                                |

## Knee age model

Image processing and model training, which produced the age predictions, are in [`code/images/`](code/images/README.md) and [`code/model/`](code/model/README.md); they are not needed to reproduce the analyses.

# 
