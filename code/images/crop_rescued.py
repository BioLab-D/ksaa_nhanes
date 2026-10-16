"""Crop ROIs of the films rescued by rotation (boxes_rot.csv).

Box coordinates refer to the rotated image, so each TIFF is first rotated with the same PIL rotate() call as in
rotate_rescue.py (angle in the rot column).
"""
import os, csv, numpy as np
from multiprocessing import Pool
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
SRC, TIF, OUT = "/home/data/knee_yolo", "/home/data/knee_xray", "/home/data/knee_roi"
MARGIN, PX = 1.15, 512

boxes = {}
for r in csv.DictReader(open(f"{SRC}/boxes_rot.csv")):
    boxes.setdefault(r["stem"], []).append(r)

def one(stem):
    rs = boxes[stem]
    try:
        im = Image.open(f"{TIF}/{stem}.tiff")
        a = np.asarray(im.rotate(int(rs[0]["rot"]), expand=True))
        H, W = a.shape
        n = 0
        for r in rs:
            dst = f"{OUT}/{stem}_{r['side']}.png"
            if os.path.exists(dst): continue
            cx, cy = float(r["cx"])*W, float(r["cy"])*H
            half = float(r["w"])*W*MARGIN/2
            x0,x1,y0,y1 = int(cx-half), int(cx+half), int(cy-half), int(cy+half)
            pad = [max(0,-y0), max(0,y1-H), max(0,-x0), max(0,x1-W)]
            roi = a[max(0,y0):min(H,y1), max(0,x0):min(W,x1)]
            if any(pad): roi = np.pad(roi, ((pad[0],pad[1]),(pad[2],pad[3])), mode="edge")
            Image.fromarray(roi.astype(np.uint16), mode="I;16").resize((PX,PX), Image.LANCZOS).save(dst)
            n += 1
        return n
    except Exception as e:
        return f"err {stem} {type(e).__name__}"

if __name__ == "__main__":
    stems = sorted(boxes)
    print(f"{len(stems)} rotated films to crop", flush=True)
    with Pool(10) as p: res = p.map(one, stems, chunksize=4)
    ok = sum(x for x in res if isinstance(x,int))
    err = [x for x in res if not isinstance(x,int)]
    print(f"new ROIs {ok}; failed {len(err)}")
    for e in err[:5]: print(" ", e)
    import glob
    print(f"total ROIs: {len(glob.glob(OUT+'/*.png'))}")
