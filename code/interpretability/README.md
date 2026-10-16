# Interpretability (Supplementary Fig. S2: occlusion sensitivity maps)

These scripts show what distinguishes older-appearing knees and where the model looks. They are not part of `run_all.sh` because they need the knee ROI images (`code/images/`), the five full-data checkpoints (`models/README.md`) and, for speed, a GPU. Their small outputs are shipped in `data/interpretability/`.

| Step | Script | Input → output (in `results/analysis/saliency/`) |
|---|---|---|
| 0 | `saliency_input.R` (R, after `run_all.sh`) | `results/analysis/analysis_data.rds` → `ksaa_for_saliency.csv` (SEQN, age, sex, out-of-fold KSAA) |
| 1 | `saliency_images.py` | ROI images → `img_all.npy` (all 4,875 participants, preprocessed as for training, 112 × 112) and `heldout_seqn.txt` (the 487 participants of the full-data models' shared internal validation set) |
| 2 | `saliency.py --heldout --bs 8` | checkpoints + ROI images → `cam_heldout.npy`, `img_heldout.npy`, `saliency_pred_heldout.csv` (Grad-CAM of predicted age, ensemble-weighted across the five backbones) |
| 3 | `saliency.py --heldout --occlude --bs 4` | → `occlusion_heldout.csv` (ensemble prediction with each of six anatomical regions blurred) |
| 4 | `saliency_figure.py` | → `saliency_summary.txt`, `saliency_regions.csv`, `results/figures/gradcam_examples.{png,pdf}` (Grad-CAM version of the example figure; not in the final manuscript) |
| 5 | `occlusion_map.py --seqn 1587,41357,32486,10939,45540,12181 --side R --tag _examples` | checkpoints + ROI images → `occlusion_map_examples.npz`: for the six example right knees, a 64 × 64-pixel patch (at 448 × 448) is blurred at 625 positions (stride 16) with the left knee unchanged; each pixel holds the mean absolute change in ensemble-predicted age over the patches covering it (about 15 s per knee on an RTX 3090) |
| 6 | `occlusion_figure.py` | → `results/figures/FigureS2.{png,pdf}` (Supplementary Fig. S2: same six knees and layout as step 4, maps lightly smoothed (σ = 8 pixels), one color scale for all panels) |
| 7 | `corner_occlusion.py --bs 1` | → `occlusion_corners_heldout.csv`, `corner_occlusion_summary.txt`: in the 487 participants, blurring each image corner, both upper corners, all four corners, or two equal-sized mid-edge control regions (both knees) changes the MAE by at most 0.14 years (both upper corners +0.09, control +0.10) |

- Paths: ROI images are read from `data/knee_roi/` (override with `KNEE_ROI=/path/to/knee_roi`) and checkpoints from `models/` (override with `KSAA_MODELS`). The model code is imported from `code/model/`.
- Only the internal validation set is used for saliency and occlusion: the full-data models were trained on everyone else and reproduce training-set ages closely (in-sample MAE about 1 year), whereas on the internal validation set their error (MAE 4.36 years) matches the cross-fitted models.
- Steps 2–3 took about 1 and 2 minutes on one NVIDIA RTX 3090 (PyTorch 2.5.1); on a laptop GPU (Apple MPS) step 2 takes about 40 minutes. CUDA and MPS results agree to within 0.002 years of predicted age.
