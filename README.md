# Knee skeletal age acceleration (KSAA) and mortality in NHANES III

Code and derived data for:

> Sheng C, Zhang Y, Li X, Chen K. *Deep learning–estimated skeletal age acceleration from knee radiographs and long-term mortality in older US adults: a cohort study.* [Journal, year, DOI — to be added]

An ensemble of five deep learning backbones estimates age from the bilateral knee radiographs of 4,875 NHANES III participants aged ≥60 years (1988–1994; ages ≥90 top-coded at 90). **Knee skeletal age acceleration (KSAA)** is the residual of out-of-fold predicted knee age regressed on chronological age, with the regression estimated in the training folds. KSAA is related to all-cause and cause-specific mortality through 2019, to Phenotypic Age acceleration and established skeletal measures, and to 482 phenotypes in a phenome-wide association study (PheWAS).

## Repository layout

```
run_all.sh            reproduces every table, figure and in-text number from data/ (about 15 min)
code/
  analysis/           R (and two Python) scripts: KSAA construction and all statistical analyses
    optional/         scripts that built shipped inputs from raw CDC files (variable dictionary, K–L grade per knee)
  images/             radiograph download, YOLO11 knee detection, ROI cropping, artifact scan   (GPU server)
  model/              knee age model adapted from Mammo-AGE; five-fold cross-fitting and full-data training (GPU server)
  interpretability/   occlusion sensitivity maps for Supplementary Fig. S2 (needs ROI images, full-data checkpoints, GPU)
data/                 inputs for code/analysis, all derived from public NHANES III files (see data/README.md)
  derived/            analysis-ready data and per-participant predictions exactly as used for the manuscript (data/derived/README.md)
models/               checksums of the full-data model checkpoints (weights distributed separately)
environment/          R and Python package versions
```

The radiographs are not redistributed. They are public at <https://ftp.cdc.gov/pub/NHANES/XRays/Nhanes3/> (4,908 knee films, about 171 GB).

## Reproducing the analyses

There are three levels:

| Level | Starts from | Needs | Reproduces |
|---|---|---|---|
| 1. Check | `data/derived/`: analysis data set, per-participant predictions, phenome-wide scan inputs and results | R | any estimate directly from the analysis-ready data; `Rscript code/analysis/check_phewas_refit.R` refits phenome-wide scan phenotypes from these files alone |
| 2. Full statistical pipeline | NHANES III files in `data/nhanes3/` and the backbones' out-of-fold predictions in `data/predictions/` | R, Python; about 15 min | every table, figure and in-text number (`bash run_all.sh`); the last step (`verify_derived`) compares the regenerated results with `data/derived/` and stops with an error if any value differs |
| 3. Images and models | CDC radiographs (171 GB) and a GPU | see `code/images/`, `code/model/` | detection, ROI cropping and model training; full-data checkpoints are distributed separately (`models/README.md`) |

### Full pipeline

Requirements: R ≥ 4.3 with the packages in `environment/R_packages.txt` (analyses used R 4.5.2); for Figure 4, a Python ≤ 3.10 environment with CLARITE 2.3.6, Pillow and PyMuPDF.

```bash
bash run_all.sh                                     # all steps; outputs go to results/
CLARITE_PY=/path/to/clarite/python bash run_all.sh  # if CLARITE is installed in a separate environment
bash run_all.sh manuscript_tables                           # re-run one step (earlier steps must have run)
```

Tables are written as HTML fragments to `results/analysis/tables/` or as CSV files to `results/analysis/weighted/`, figures to `results/figures/`, and console logs to `results/logs/`.

Primary analyses are survey-weighted (NHANES III examination weights, strata and primary sampling units; phase 2 weights for analyses of the radiographic K–L reading). `weighted_common.R` builds the design on all 17,705 adult examinees (table `nhanes3_design_all` in `data/nhanes3/nhanes3_knee_cohort.duckdb`) and analyzes the cohort as a domain. Bootstrap CIs use Rao–Wu rescaling bootstrap replicate weights (`weighted_boot.R`; 500 replicates for standardized risks, 200 for discrimination). The unweighted and multiple-imputation estimates in Table S5 are computed in `weighted.R`.

### Where each result comes from

Outputs used in the manuscript are named after the manuscript item (for example `figures/Figure2.pdf`, `analysis/weighted/tableS5.csv`, `analysis/tables/table1.html`); intermediate results have descriptive names.

| Manuscript item | Output (under `results/`) | Script(s) |
|---|---|---|
| KSAA | `ksaa_oof.csv` | `ksaa_ensemble.R` |
| Table 1 | `analysis/tables/table1.html` | `manuscript_tables.R` |
| Table 2 | `analysis/weighted/table2.csv` | `weighted.R` |
| Figure 1 | `figures/Figure1/Figure1.{pdf,png}` | `weighted.R`, `weighted_boot.R`, `phewas.R` → `fig1.R` |
| Figure 2 | `figures/Figure2.{pdf,png}`; `analysis/weighted/figure2_{std_curves,absrisk,ptrend}.csv` | `weighted_boot.R` → `fig2_mortality.R` |
| Figure 3 | `figures/Figure3.{pdf,png}`; `analysis/weighted/figure3a.csv`, `figure3_cells.csv`, `figure3b_cindex.csv` | `weighted.R`, `weighted_boot.R` → `fig3_incremental.R` |
| Figure 4 | `figures/Figure4.{pdf,png}` | `phewas.R` → `fig4_phewas.R` (panel B), `fig4_phewas_clarite.py` (panel A and assembly) |
| Table S1 | `analysis/tables/tableS1.html` | `manuscript_tables.R`, `tableS1_performance.R` |
| Table S2 | covariate definitions (text table; implemented in `build_analysis_data.R`) | `build_analysis_data.R` |
| Table S3 | `analysis/tables/tableS3.html` (means and percentages unweighted; P values from design-based Rao–Scott and Kruskal–Wallis tests) | `quartile_characteristics.R` → `manuscript_tables.R` |
| Table S4 | `analysis/weighted/tableS4.csv` | `weighted.R` |
| Table S5 | `analysis/weighted/tableS5.csv` (also contains rows not shown in the manuscript) | `weighted.R` |
| Table S6 | `analysis/weighted/tableS6.csv` | `weighted_boot.R` |
| Table S7 | `analysis/weighted/tableS7.csv` | `weighted.R` |
| Table S8 | `analysis/tables/tableS8.html` | `phewas_clean.R` (`EWAS_STAGE=audit`), `phewas_audit.R`, `phewas.R` → `phewas_table.R` |
| Table S9 | `analysis/tables/tableS9.html`; all results in `analysis/phewas/phewas_results.csv` | as for Table S8 |
| Table S10 | `analysis/weighted/tableS10.csv` (also contains a combined-fracture row and a "retained" column not shown in the manuscript) | `weighted.R` |
| Fig. S1 | `figures/FigureS1.{pdf,png}` | `figS1_biascorr.R` |
| Fig. S2 | `figures/FigureS2.{pdf,png}`, occlusion sensitivity maps (run separately; see `code/interpretability/README.md`) | `code/interpretability/occlusion_map.py` → `occlusion_figure.py` |
| Fig. S3 | `figures/FigureS3.{pdf,png}`; `analysis/weighted/figureS3.csv` | `weighted.R` → `fig2_mortality.R` |
| Partial correlations of KSAA with the K–L grade and femoral neck BMD (Section 3.4) | `analysis/kl_dxa_partial_cor.csv` (columns `r`, `r_lo`, `r_hi`) | `kl_bmd_partial_cor.R` |
| In-text numbers (follow-up, person-years, etc.) | `logs/manuscript_numbers.txt` and the CSV files in `analysis/` and `analysis/weighted/` | `manuscript_numbers.R` |

The Figure 2 P values for trend use all 4,875 participants (age and sex adjusted); the Table 2 P values use the 4,347 participants with complete Model 3 covariates, so the two differ.

Outputs in `results/` are manuscript items or inputs to them, plus console logs (`logs/`) and working copies of inputs (`work/`, `predictions/`). Analyses that were reported only in earlier drafts (unweighted versions of the main analyses, comparisons with other markers, correlates, bone and mineral markers, supplementary figures and tables that were later dropped) have been removed from this repository.

## Upstream steps (not needed to reproduce the analyses)

1. **Images** (`code/images/`): download the CDC TIFFs, convert to 8-bit PNG for detection, fine-tune YOLO11s on bootstrapped labels, detect both knees, crop a square ROI (1.15 × the bounding-box width, 16-bit), recover failed detections by rotation search or relaxed thresholds, and scan ROIs for bright artifacts. The outputs used downstream are in `data/images/`.
2. **Model** (`code/model/`): `knee_v2.py` trains the adapted Mammo-AGE model (shared ImageNet-pretrained backbone for both knees, four-block transformer, probabilistic ordinal embedding heads over 100 one-year age classes, POE and mean–variance losses). `train_5fold.sh` (section 1) runs the five pre-specified backbones (ResNet-18, ResNet-50, EfficientNet-B0, ConvNeXt-Tiny, DenseNet-121) under participant-level five-fold cross-fitting (split seed 20260831); later sections of the script are ablations not reported in the manuscript. The out-of-fold predictions are in `data/predictions/`. `train_full.sh` trains each backbone once on all participants for use on new radiographs.
3. **Interpretability** (`code/interpretability/`): occlusion sensitivity maps (Fig. S2) and Grad-CAM maps of the full-data models in their internal validation set (487 participants not used for training), plus region and corner occlusion tests. See `code/interpretability/README.md`.

## Data sources and terms of use

All inputs are derived from public-use files of the National Center for Health Statistics (NCHS): the NHANES III (1988–1994) household adult, examination and laboratory files; the NHANES III knee radiographs; and the NCHS 2019 public-use Linked Mortality File. Users must follow the NCHS data-use restrictions, including that no attempt be made to identify any participant. `data/README.md` describes every shipped file.

## Licenses

- Code written for this study: [license to be added].
- `code/model/mammo_age_model.py`, `transformer.py`, `probordiloss.py`, `mean_variance_loss.py` and `metrics.py` are from the official Mammo-AGE implementation (Apache License 2.0, `code/model/LICENSE_Mammo-AGE`). `mammo_age_model.py` and `transformer.py` were modified: the number of views is set by the environment variable `N_VIEWS` (2 knees instead of 4 mammographic views), the global dropout rate by `DROP`, and the import path of the transformer module was changed; the other three files are unchanged.
- Knee detection uses Ultralytics YOLO11 (AGPL-3.0) as an installed dependency; no Ultralytics code or trained detector weights are included.
- NHANES data are in the public domain.

## Notes

- Scripts are the versions that produced the manuscript, with paths made relative to this repository, version labels removed from file names, comments condensed, and a few guards for legacy files that are not shipped. Many code comments are in Chinese.
- Follow-up is measured from the mobile examination center visit (`PERMTH_EXM`) to death or 31 December 2019.
- Phenotypic Age acceleration is winsorized at the 1st and 99th percentiles in the analyses.

## Citation

If you use this code or the derived data, please cite the article above. [Citation to be added on publication]

## Contact

Xiangchun Li (lixiangchun@tmu.edu.cn) and Kexin Chen (chenkexin@tmu.edu.cn), Tianjin Medical University Cancer Institute and Hospital.
