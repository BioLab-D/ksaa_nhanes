"""Crop ROIs of the films rescued with relaxed detection thresholds (boxes_lowconf.csv).

boxes_lowconf.csv has one row per film with columns for each knee, and coordinates refer to the unrotated TIFF.
Cropping follows crop_rois.py (square box 1.15 x bone width, 512 px, 16-bit), so these 107 films have the same ROI scale
as the rest of the sample.
"""
import os, csv, numpy as np
from multiprocessing import Pool
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
SRC, TIF, OUT = "/home/data/knee_yolo", "/home/data/knee_xray", "/home/data/knee_roi"
MARGIN, PX = 1.15, 512

rows = list(csv.DictReader(open(f"{SRC}/boxes_lowconf.csv")))

def one(r):
    stem = r["stem"]
    try:
        a = np.asarray(Image.open(f"{TIF}/{stem}.tiff")); H, W = a.shape
        n = 0
        for side in ("R", "L"):
            dst = f"{OUT}/{stem}_{side}.png"
            if os.path.exists(dst): continue
            cx, cy = float(r[f"cx{side}"])*W, float(r[f"cy{side}"])*H
            half = float(r[f"w{side}"])*W*MARGIN/2
            x0,x1,y0,y1 = int(cx-half), int(cx+half), int(cy-half), int(cy+half)
            pad = [max(0,-y0), max(0,y1-H), max(0,-x0), max(0,x1-W)]
            roi = a[max(0,y0):min(H,y1), max(0,x0):min(W,x1)]
            if any(pad): roi = np.pad(roi, ((pad[0],pad[1]),(pad[2],pad[3])), mode="edge")
            Image.fromarray(roi.astype(np.uint16), mode="I;16").resize((PX,PX), Image.LANCZOS).save(dst)
            n += 1
        return n
    except Exception as e:
        return f"err {stem} {type(e).__name__}: {e}"

if __name__ == "__main__":
    print(f"{len(rows)} low-confidence films to crop", flush=True)
    with Pool(10) as p: res = p.map(one, rows, chunksize=4)
    ok  = sum(x for x in res if isinstance(x, int))
    err = [x for x in res if not isinstance(x, int)]
    print(f"new ROIs {ok}; failed {len(err)}", flush=True)
    for e in err[:10]: print(" ", e)
