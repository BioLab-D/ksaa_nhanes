# Knee age model

`knee_v2.py` adapts Mammo-AGE (Wang et al., Nat Commun 2025) by treating a participant's two knees as the two views. Each knee passes through a shared ImageNet-initialised backbone. An instance-bag transformer (`--nblock 4`) fuses the two knees. Knee-level probabilistic ordinal embedding heads cover 100 one-year age classes and are trained with a probabilistic ordinal loss and a mean–variance loss. The prediction is the mean of the two knee-level expected ages.

Configuration used for the manuscript (see `train_5fold.sh`, section 1):

```
python knee_v2.py --tag cv_<backbone> --arch <resnet18|resnet50|efficientnet_b0|convnext_tiny|densenet121> --img 448 \
  --bs 16 --seed 20260831 --folds 0,1,2,3,4 --nblock 4 \
  --aug knee --norm pct --sched cosine --wd 1e-4 --warmup 5 --epochs 80 --min_ep 60 --patience 15
# environment: N_VIEWS=2 DROP=0
```

- Participant-level outer folds use `--split_seed 20260831` (default). Within each outer fold, 10% of the training participants form the internal validation set used for early stopping and the ensemble weights. The held-out fold is used only for prediction.
- Knee-specific recipe: AdamW with cosine decay (learning rate 1 × 10⁻⁴, weight decay 1 × 10⁻⁴, five-epoch warm-up), up to 80 epochs, and early stopping on internal-validation MAE (patience 15, not before epoch 60). Augmentations and normalisation are listed in the docstring of `knee_v2.py`.
- Paths are those of the training server: the shell scripts run `knee_v2.py` from `/tmp/mammo_age_official` (a copy of this directory) and write to `/home/data/knee_yolo/v2`; edit them for your system.
- Inputs: the ROI directory (`--roi_dir`, default `/home/data/knee_roi`) and a label file `/tmp/surv_labels.csv` with columns `SEQN` and `age` (produced by `code/analysis/01_survival_labels.R`).
- Outputs: `{tag}_oof.csv` (out-of-fold predictions) and `{tag}_folds.csv` (per-fold internal-validation MAE), as shipped in `data/predictions/`. Held-out predictions were written to a separate directory and were not examined until all backbones had been trained on all folds.
- `train_full.sh`: `--full --save_dir <dir>` trains each backbone once on all participants and saves the checkpoint (see `models/README.md`).
- Hardware: one NVIDIA RTX 3090 (24 GB). Mean training time per outer fold ranged from 111 min (ResNet-18) to 424 min (ConvNeXt-Tiny), about 100 GPU-hours for the five backbones × five folds.

The ensemble, bias correction and KSAA are computed in R (`code/analysis/02_ksaa_ensemble.R`).
