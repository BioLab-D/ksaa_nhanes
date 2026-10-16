"""Train the knee ROI detector (2 classes: participant's right knee / left knee).

Labels come from heuristic bootstrapping, not a manual gold standard - 23 of 24 checked boxes were correct (~96%).
Remaining noise is handled in two ways:
  1. automatic exclusion by suspicion score before training (weak joint-localization dip, or left-right joint-row difference >250)
  2. the detector averages over label noise and learns anatomy rather than memorizing films

Key setting fliplr=0.0: the knees are mirror images, so horizontal-flip augmentation would make the two classes
indistinguishable in pixels and break the laterality head (same reason flips are off in eus_det/train_five.py).
"""
import os, csv, random, shutil, sys
SRC = "/home/data/knee_yolo"
DS  = "/home/data/knee_yolo/ds"

rows = [r for r in csv.DictReader(open(f"{SRC}/meta.csv")) if r["ok"] == "1"]
clean = [r for r in rows if float(r["susp"]) == 0.0]
print(f"total {len(rows)}  clean {len(clean)}  excluded {len(rows)-len(clean)}")

random.seed(0); random.shuffle(clean)
n_val = max(20, int(len(clean)*0.2))
split = {"val": clean[:n_val], "train": clean[n_val:]}

for sp in ("train", "val"):
    for kind in ("images", "labels"):
        d = f"{DS}/{kind}/{sp}"
        shutil.rmtree(d, ignore_errors=True); os.makedirs(d, exist_ok=True)
    for r in split[sp]:
        os.link(f"{SRC}/images/{r['stem']}.png", f"{DS}/images/{sp}/{r['stem']}.png")
        os.link(f"{SRC}/labels/{r['stem']}.txt", f"{DS}/labels/{sp}/{r['stem']}.txt")
    print(f"  {sp}: {len(split[sp])}")

with open(f"{DS}/data.yaml", "w") as f:
    f.write(f"path: {DS}\ntrain: images/train\nval: images/val\n"
            "names:\n  0: knee_R\n  1: knee_L\n")

if "--train" in sys.argv:
    from ultralytics import YOLO
    m = YOLO("yolo11s.pt")
    m.train(data=f"{DS}/data.yaml", epochs=80, imgsz=1024, batch=8,
            device=0, workers=8, patience=25,
            project="/home/data/knee_yolo/runs", name="det_v1", exist_ok=True,
            fliplr=0.0,      # knees are mirror-symmetric; never flip horizontally
            flipud=0.0,
            degrees=4.0, translate=0.08, scale=0.25, shear=0.0,
            mosaic=0.5, close_mosaic=15, hsv_h=0.0, hsv_s=0.0, hsv_v=0.3)
