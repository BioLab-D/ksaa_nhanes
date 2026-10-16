# Trained models

The manuscript's KSAA values come from the five-fold cross-fitted models (`code/model/train_5fold.sh`); their out-of-fold predictions are in `data/predictions/`. Fold-level weights were not retained.

For application to new knee radiographs, each of the five pre-specified backbones was also trained once on all 4,875 participants (`code/model/train_full.sh`; 10% internal validation set for early stopping only). The checkpoints (about 586 MB in total; two exceed 100 MB) are distributed separately:

- Location: [Zenodo DOI — to be added]
- Files and SHA-256 checksums (check with `shasum -a 256 m5_full_*_full.pt`):

```
f5d8e1429f731ba14899ddfcef2bd23c491452dc6b49fa7f90a0d5c3ab88f3d3  m5_full_convnexttiny_full.pt
760b8fc5701546d6be19237b10d7635dac9eea3fac009be5c119a902ac356f2d  m5_full_densenet121_full.pt
528928477c77d2839c3a6b918db118381173f4bd7c4ba86d5303c9f73245bc57  m5_full_efficientnetb0_full.pt
e12871700686af0ced26f9d0b379a7c0ba72b33bce13fb2da3ad8d7f39b7a583  m5_full_resnet18_full.pt
5c3e8a1b961b9a358f6cdad3dca3a6320326d5270e96dfe2dba1c73ab169c371  m5_full_resnet50_full.pt
```

Each checkpoint (`m5_full_<backbone>_full.pt`) is a dict with `state_dict`, `arch`, the intensity-normalization parameters (`par`), all training arguments (`args`), `views`, `min_age` and the internal-validation participants (`val_seqn`). Rebuild the network with `Mammo_AGE` from `code/model/mammo_age_model.py` (environment variable `N_VIEWS=2`) using `arch` and the `num_output`, `nblock` and `hidden` entries of `args`, then load `state_dict` (see `code/interpretability/saliency.py`); the predicted age of one checkpoint is the mean of the two knee-level expected ages. If you have the files, place them in this folder: `.gitignore` excludes `models/*.pt` from the repository.

## Applying the models to new radiographs

As stated in each checkpoint's `note`, combine the five backbones with the ensemble weights and bias correction calibrated on the five-fold out-of-fold predictions, given in `data/derived/new_data_calibration.csv`. Do not re-estimate them from the full-data models' predictions for the 4,875 training participants, which are in-sample and would understate the age bias.

```
pred_raw = sum over the five checkpoints of ensemble_weight × predicted age        # weights: means of the fold-specific weights, sum to 1
KSAA     = pred_raw − (expected_pred_intercept + expected_pred_slope × age)         # 33.555 + 0.5271 × age, fitted on all out-of-fold predictions
KSAA per SD = KSAA / ksaa_sd                                                        # 3.945 years, the SD in the study cohort
```

The fold-specific weights and bias-correction coefficients used for the study's own (out-of-fold) KSAA are in `data/derived/ensemble_weights_by_fold.csv` and `bias_correction_by_fold.csv`. Predictions for radiographs from other sources (weight-bearing views, digital radiography, adults younger than 60 years) have not been validated.
