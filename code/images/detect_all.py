"""Run the knee detector on all 4,908 films and write box coordinates and a quality category.

Boxes are stored as fractions of the film size, so they map directly onto the original TIFF.
A film is "ok" when exactly one right and one left knee are detected side by side; all others go to the review list.
"""
import os, glob, csv, numpy as np
from ultralytics import YOLO
SRC = "/home/data/knee_yolo"

paths = sorted(glob.glob(f"{SRC}/images/*.png"))
print(f"{len(paths)} films to detect", flush=True)
m = YOLO(f"{SRC}/runs/det_v1/weights/best.pt")

fb = open(f"{SRC}/boxes.csv","w",newline=""); wb = csv.writer(fb)
wb.writerow(["stem","side","cx","cy","w","h","conf"])          # coordinates as fractions of film size
fq = open(f"{SRC}/qc.csv","w",newline=""); wq = csv.writer(fq)
wq.writerow(["stem","status","n_R","n_L","dy_norm"])

cat = {}
for i in range(0, len(paths), 8):
    for p, rr in zip(paths[i:i+8],
                     m.predict(paths[i:i+8], imgsz=1024, conf=0.25,
                               device=0, verbose=False)):
        stem = os.path.basename(p)[:-4]
        H, W = rr.orig_shape
        b = rr.boxes
        cls = b.cls.cpu().numpy().astype(int) if b is not None and len(b) else np.array([])
        xywh = b.xywh.cpu().numpy() if b is not None and len(b) else np.zeros((0,4))
        cf = b.conf.cpu().numpy() if b is not None and len(b) else np.array([])
        nR, nL = int((cls==0).sum()), int((cls==1).sum())
        dy = ""
        if nR==1 and nL==1:
            pR, pL = xywh[cls==0][0], xywh[cls==1][0]
            dy = round(float(abs(pR[1]-pL[1])/H), 4)
            st = "ok" if abs(pR[0]-pL[0]) >= abs(pR[1]-pL[1]) else "rotated"
            if st=="ok":
                for c, xy, cc in zip(cls, xywh, cf):
                    wb.writerow([stem, "R" if c==0 else "L",
                                 round(xy[0]/W,6), round(xy[1]/H,6),
                                 round(xy[2]/W,6), round(xy[3]/H,6), round(float(cc),4)])
        elif nR+nL==0: st = "none"
        elif nR>1 or nL>1: st = "dup"
        else: st = "single"
        wq.writerow([stem, st, nR, nL, dy]); cat[st] = cat.get(st,0)+1
    if (i//8) % 50 == 0: print(f"  {i}/{len(paths)}", flush=True)
fb.close(); fq.close()

N = sum(cat.values())
print(f"\n=== {N} films ===")
for k,v in sorted(cat.items(), key=lambda x:-x[1]):
    print(f"  {k:<9} {v:>5}  ({v/N:6.2%})")
print(f"\nboxes.csv written: {cat.get('ok',0)*2} boxes")
