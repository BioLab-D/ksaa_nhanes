# Derived analysis data

Analysis-ready data and per-participant predictions, exactly as used for the manuscript. They let readers check any result directly, without running the full pipeline, and they are the reference that `run_all.sh` compares its regenerated outputs against (step `verify_derived`, script `code/analysis/verify_derived.R`).

All values are derived from public-use NCHS files (NHANES III, 1988–1994; NHANES III knee radiographs; public-use Linked Mortality File, follow-up through 31 December 2019). `SEQN` links all files; rows are in the same order in every per-participant file. Age is top-coded at 90 years. Missing values are empty cells.

| File | Content | Produced by |
|---|---|---|
| `analysis_dataset.csv` | 4,875 participants × 84 columns: KSAA, outcomes, survey design, covariates, Phenotypic Age and its biomarkers, musculoskeletal variables. Read by all mortality, descriptive and subgroup analyses | `code/analysis/03_build_analysis_data.R` |
| `data_dictionary.csv` | Meaning, unit or coding, and NHANES III source variable of each column of `analysis_dataset.csv` | — |
| `ensemble_predictions.csv` | Per participant: outer fold, chronological age, `pred_raw` (weighted ensemble of the five backbones), `pred_corr` (bias-corrected), `gap_raw`, `ksaa`, `q` (KSAA quartile) | `code/analysis/02_ksaa_ensemble.R` |
| `ensemble_weights_by_fold.csv` | Weight of each backbone in each outer fold: 1 / internal-validation MAE, normalized within the fold | from `data/predictions/*_folds.csv` |
| `bias_correction_by_fold.csv` | Age-bias correction for each held-out fold: `pred_raw − age ~ age` fitted on the other four folds | from `ensemble_predictions.csv` |
| `new_data_calibration.csv` | Ensemble weights (mean of the fold-specific weights) and bias correction (`pred_raw ~ age` fitted on all out-of-fold predictions) for applying the full-data checkpoints in `models/` to new radiographs, and the SD of KSAA; see `models/README.md` | from `ensemble_predictions.csv` and `ensemble_weights_by_fold.csv` |
| `phewas/variable_audit_table.csv` | Variable audit (Table S8): health domain, type, eligible population, survey phase, valid-response rate, counts, status (main / supplementary / exposure_module / covariate / excluded) and reason, summaries of repeated measurements and recodings | `code/analysis/09_phewas_audit.R` |
| `phewas/phewas_results.csv` | Results for all 715 phenotypes (main scan 482, supplementary 59, exploratory module 174; Table S9, Figure 4) | `code/analysis/10_phewas.R` |
| `phewas/phenotype_values.csv` | The value of each of 712 phenotypes **as entered into the regression** (continuous: log-transformed if right-skewed, then standardized; binary: 0/1; ordinal: level index after merging sparse levels, see `levels` in the results). Empty where the participant was not in that analysis. The three nominal phenotypes of the exploratory module are not included | `code/analysis/10_phewas.R` |
| `phewas/survey_design.csv` | Six-year, phase-specific and morning fasting-subsample weights, PSUs and strata used by the scan (missing weights set to 0, as in the scan) | NHANES III |
| `phewas/phewas_summary.txt` | Counts per set, independent signals and the top phenotypes | `code/analysis/10_phewas.R` |

The out-of-fold predictions of each backbone are in `data/predictions/`, and the 825-variable candidate pool and the multilevel plan are in `data/phewas/`.

## Recomputing KSAA from the backbone predictions

```
pred_raw = sum over backbones of weight[fold, backbone] × pred_age[backbone]     # ensemble_weights_by_fold.csv, data/predictions/cv_<backbone>_oof.csv
ksaa     = pred_raw − age − (intercept[fold] + slope[fold] × age)                # bias_correction_by_fold.csv
```

`ksaa` equals `gap` (and `gap_corr`) in `analysis_dataset.csv`.

## Phenome-wide scan

Each phenotype was regressed on KSAA per SD (`gap / sd(gap)`, SD 3.95 years) with adjustment for age, sex, race and ethnicity, education (`edu3`) and poverty–income ratio, using `WTPFEX6`, `SDPPSU6` and `SDPSTRA6` with `nest = TRUE` (morning fasting-subsample weights for fasting assays; the phase design for single-phase items). Continuous: `svyglm`; binary: `svyglm` with quasibinomial family; ordinal: `svyolr` (refitted as any level above the lowest vs the lowest when singular); nominal: joint test from a regression of KSAA on the categories (P only). P values and 95% CIs use a t distribution with the design degrees of freedom. FDR is controlled within each `set`. `clump_lead` is the lead phenotype of the independent signal. The columns `bonf`, `meff` and `p_liji` (Bonferroni and Li–Ji corrections) are not used in the manuscript.

`code/analysis/check_phewas_refit.R` refits a set of phenotypes from these files alone (R with `survey` and `MASS`) and confirms that n, estimates and P values match `phewas_results.csv`.

## Notes

- Survey-weighted analyses in the manuscript define the design on all adult examinees and treat the cohort as a domain (`code/analysis/weighted_common.R`); `options(survey.lonely.psu = "adjust")`.
- Follow-up time `t` starts at the mobile examination center examination (`PERMTH_EXM`), in months.
- `gap_L` and `gap_R` come from the left and right knee-level heads of the bilateral model. These heads also receive information from the other knee (cross-attention and a shared bilateral branch), so they are not independent single-knee estimates.
