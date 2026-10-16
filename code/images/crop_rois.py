"""Crop knee ROIs from the original 16-bit TIFFs using the detection boxes.

ROIs are kept at 16 bits so that the grey-level window can be chosen at training time rather than fixed here.
The box is enlarged 1.15 x for anatomical context; regions beyond the film edge are edge-padded, not zero-padded,
to avoid artificial edges at the ROI border.
"""
import os, csv, numpy as np
from multiprocessing import Pool
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
SRC, TIF, OUT = "/home/data/knee_yolo", "/home/data/knee_xray", "/home/data/knee_roi"
MARGIN, PX = 1.15, 512

boxes = {}
for r in csv.DictReader(open(f"{SRC}/boxes.csv")):
    boxes.setdefault(r["stem"], []).append(r)

def one(stem):
    rs = boxes[stem]
    try:
        a = np.asarray(Image.open(f"{TIF}/{stem}.tiff"))
        H, W = a.shape
        out = []
        for r in rs:
            side = r["side"]
            dst = f"{OUT}/{stem}_{side}.png"
            if os.path.exists(dst): out.append("skip"); continue
            cx, cy = float(r["cx"])*W, float(r["cy"])*H
            half = float(r["w"])*W*MARGIN/2
            x0, x1 = int(cx-half), int(cx+half)
            y0, y1 = int(cy-half), int(cy+half)
            pad = [max(0,-y0), max(0,y1-H), max(0,-x0), max(0,x1-W)]
            roi = a[max(0,y0):min(H,y1), max(0,x0):min(W,x1)]
            if any(pad):
                roi = np.pad(roi, ((pad[0],pad[1]),(pad[2],pad[3])), mode="edge")
            Image.fromarray(roi.astype(np.uint16), mode="I;16") \
                 .resize((PX,PX), Image.LANCZOS).save(dst)
            out.append("ok")
        return out
    except Exception as e:
        return [f"err {stem} {type(e).__name__}: {e}"]

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    stems = sorted(boxes)
    print(f"{len(stems)} films / {sum(len(v) for v in boxes.values())} knees to crop", flush=True)
    with Pool(12) as pool:
        res = pool.map(one, stems, chunksize=8)
    flat = [x for sub in res for x in sub]
    from collections import Counter
    print(dict(Counter(x if x in ("ok","skip") else "err" for x in flat)))
    for x in flat:
        if x not in ("ok","skip"): print(" ", x)
