# Data

Every file here is derived from public-use NCHS files. `SEQN` is the NHANES III respondent sequence number and links all files. Variable names are the original NHANES III names unless stated; missing-value codes (8, 88, 888, …) are left as in the source files and are removed variable by variable in the scripts.

Sources:
- NHANES III household adult (`adult.dat`), examination (`exam.dat`) and laboratory (`lab.dat`) files, 1988–1994: https://wwwn.cdc.gov/nchs/nhanes/nhanes3/
- NHANES III knee radiograph reading file (`xr.dat`, phase 2) and the knee radiographs (released 2025): https://ftp.cdc.gov/pub/NHANES/XRays/Nhanes3/
- NCHS public-use Linked Mortality File, NHANES III, follow-up through 31 December 2019: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

## nhanes3/

`nhanes3_knee_cohort.duckdb` (DuckDB 1.x) — the complete NHANES III adult (1,239 columns), examination (2,369), laboratory (357) and linked-mortality (14) files, restricted to the 4,908 participants with released knee radiographs (table `cohort_seqn`), plus `nhanes3_value_labels`. Read by `01_survival_labels.R` (survival labels), `03_build_analysis_data.R` (covariates, Phenotypic Age biomarkers), `08_phewas_clean.R` and `10_phewas.R` (phenome-wide scan).

## predictions/

Output of `code/model/train_5fold.sh` (section 1), one pair of files per pre-specified backbone (`resnet18`, `resnet50`, `efficientnetb0`, `convnexttiny`, `densenet121`):

- `cv_<backbone>_oof.csv` — out-of-fold predictions for all 4,875 participants: `SEQN`, `fold` (outer fold 0–4), `age` (chronological age, years), `pred_age` (predicted age used in the analyses: mean of the two knee-level expected ages), `pred_global`, `pred_local`, `pred_L`, `pred_R` (other read-outs of the same model, not used).
- `cv_<backbone>_folds.csv` — per outer fold: `inner_val_mae` (mean absolute error in the internal validation set drawn from the training folds; used for the fold-wise ensemble weights), `best_ep`, `stop_ep`, `minutes`.

## prepared/

Variable extracts from the NHANES III files (all participants of the respective file unless stated), created during the project with the NHANES III SAS layouts:

| File | Content |
|---|---|
| `demo_all.csv` | age, sex, race and ethnicity (adult file) |
| `dxa_all.csv` | proximal femur DXA (`BDP*`: BMD, BMC and area for femoral neck, trochanter, intertrochanter, Ward's triangle; total) |
| `fracture.csv` | self-reported fractures (`HAG5–HAG15`) |
| `knee_pain.csv` | knee pain items (`HAG22–HAG34`) and selected health-history items (`HAC1A`, `HAC1D`, `HAD1`, `HAF10`, `HAR1`, `HAR3`) |
| `physfun.csv` | physical performance (8-foot walk, tandem stand, `PFP*`) and functional difficulty (`HAH*`) |
| `gait.csv` | 8-foot walk times (`PFPTWLKA`, `PFPTWLKB`) and related examination items |
| `validate_vars.csv` | comorbidity and smoking items (`HAD1`, `HAC1A`, `HAC1D`, `HAF10`, `HAR1`) and other health-history variables; read by `02_ksaa_ensemble.R` for the covariates of its Cox models |
| `kl_per_knee.csv` | Kellgren–Lawrence grade per knee from the radiograph reading (`xr.dat`, phase 2): `seqn`, `side`, `kl`, `oa` (K–L ≥ 2) |
| `xray_components.csv` | all scored features of the radiograph reading (`XRP*`: K–L grade, osteophyte and sclerosis scores by compartment, chondrocalcinosis, joint replacement), phase 2 |

## derived/

Analysis-ready data and per-participant predictions exactly as used for the manuscript, and the reference for the final check of `run_all.sh`; see `derived/README.md`.

## images/

| File | Content |
|---|---|
| `cdc_knee_seqn.txt` | SEQN of the 4,908 participants with released knee radiographs |
| `boxes_rot.csv` | detections recovered by rotation search (163 films): `stem` (file stem, `<SEQN>K1`), rotation, side, box centre/width/height, confidence |
| `boxes_lowconf.csv` | detections recovered with relaxed thresholds (107 films) |
| `metal_scan.csv` | per-ROI bright-artefact scan: `roi`, `frac_top` (fraction of pixels near the maximum), `p999_ratio`; used to flag films with bright artifacts (exploratory; not reported in the manuscript) |
| `example_roi/10018K1_{L,R}.png` | the two 16-bit knee ROIs of one participant, shown in Figure 1A |

## interpretability/

Outputs of `code/interpretability/`. Supplementary Fig. S2 shows the occlusion sensitivity maps in `occlusion_map_examples.npz`; the Grad-CAM, region and corner occlusion summaries are additional and not reported in the manuscript. The saliency and occlusion files cover the 487 participants of the full-data models' shared internal validation set, who were not used to train those models.

| File | Content |
|---|---|
| `heldout_seqn.txt` | SEQN of the 487 internal-validation participants |
| `saliency_pred_heldout.csv` | `SEQN`, `age`, `gap` (out-of-fold KSAA), `pred_oof` (uncorrected out-of-fold ensemble prediction), predicted age of each full-data backbone and their weighted ensemble (`pred_full_ens`) |
| `occlusion_heldout.csv` | the same predictions with no region blurred (`region` = `None`) and with each of six regions blurred in both knees |
| `occlusion_map_examples.npz` | occlusion sensitivity maps of the six example right knees shown in Fig. S2 (`seqn`, `abs` = mean absolute change in predicted age, `signed`, `base` = unoccluded prediction; 112 × 112) |
| `occlusion_corners_heldout.csv`, `corner_occlusion_summary.txt` | corner occlusion test: predictions with each image corner, both upper corners, all four corners, or two mid-edge control regions blurred in both knees, and the change in MAE with 95% bootstrap CIs |
| `saliency_regions.csv` | per region: share of image area and of Grad-CAM saliency, enrichment with 95% bootstrap CI, occlusion change in MAE and in predicted age with 95% bootstrap CIs, and the mean of the age and KSAA difference maps |
| `saliency_summary.txt` | summary of the saliency and occlusion analyses |

## phewas/

Dictionary of the 3,962 NHANES III adult, examination and laboratory variables, built by `code/analysis/optional/phewas_dict.R` from the CDC SAS layout files: `phewas_vars.csv` (position, label, section), `phewas_codes.csv` (code ranges and value labels, including missing-value codes), `phewas_dict_check.csv` (checks of the codes against the observed values). These drive the variable audit (Supplementary Table S8) and the phenome-wide scan.

Two further files fix the inputs of the scan as they were used for the manuscript:

- `phewas_candidates.csv` — the 825-variable candidate pool (`analysed = TRUE`) produced by a preliminary scan with an earlier version of `code/analysis/08_phewas_clean.R`, to which the pre-specified audit protocol (`09_phewas_audit.R`) was then applied. After the pool was fixed, `08_phewas_clean.R` was changed to treat fundus-photograph "Can't grade" codes as missing at the candidate stage; with that change six fundus grading items (`FPP1110`–`FPP1160`) would drop out of the pool, so the pool is shipped rather than regenerated.
- `phewas_multilevel_plan.csv` — the classification of multi-category phenotypes (ordinal or nominal) used by `10_phewas.R`, produced by `optional/phewas_multilevel_plan.R` from the preliminary scan.
