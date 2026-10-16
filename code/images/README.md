# Radiograph processing

These scripts ran on a Linux GPU server with the data under `/home/data/` (paths are hard-coded at the top of each script, and `pipeline.sh` calls the Python scripts from `/tmp`; edit them for your system). They are not needed to reproduce the statistical analyses, whose image-derived inputs are shipped in `data/images/` and `data/predictions/`.

Order of execution:

| Step | Script | Input → output |
|---|---|---|
| 1 | `knee_xray_download.sh` | CDC FTP (`https://ftp.cdc.gov/pub/NHANES/XRays/Nhanes3/`) → 4,908 knee TIFFs (16-bit, about 5000 × 4200 px, 171 GB). The script resumes interrupted downloads and checks each file against the byte size CDC lists; also check the TIFF magic number, because an intercepting portal on some networks returns HTML pages of the expected size. |
| 2 | `convert_all.py` | TIFF → 1024-px 8-bit PNG (detector input only) |
| 3 | `make_yolo_ds.py`, `train_knee.py` | heuristic bootstrap of knee bounding boxes and fine-tuning of YOLO11s (two classes: participant's right and left knee; horizontal flip disabled) |
| 4 | `pipeline.sh` → `detect_all.py`, `crop_rois.py` | detection on all films; square ROI centred on each knee with side 1.15 × detected bone width, cropped from the original 16-bit TIFF and saved as 512-px 16-bit PNG (`<SEQN>K1_{L,R}.png`) |
| 5 | `rotate_rescue.py`, `crop_rescued.py` | films without a valid detection: rotation search (±90°) and cropping (163 films; `data/images/boxes_rot.csv`) |
| 6 | `crop_lowconf.py` | remaining failures accepted with relaxed detection thresholds (107 films; `data/images/boxes_lowconf.csv`) |
| 7 | `metal_scan.py` | bright-artefact scan of all ROIs (`data/images/metal_scan.csv`) |

33 films (0.7%) could not be processed; 4,875 participants have both knee ROIs. The model mirrors left knees and normalises intensities at training time (`code/model/knee_v2.py`), so the saved ROIs are unmirrored and unnormalised.
