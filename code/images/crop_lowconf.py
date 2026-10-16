"""裁低置信度救援片(boxes_lowconf.csv)。

与 boxes_rot.csv 不同,这批是**宽表**:一行一张片,左右膝各自成列,
且没有旋转 —— 坐标直接相对原始 TIFF。裁法与 crop_rois.py 完全一致
(骨宽 1.15 倍方框、512px、16-bit),否则这 107 张的 ROI 尺度会和
主样本的 9,536 个对不齐,合进去等于引入新的批次变量。
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
    print(f"待裁 {len(rows)} 张低置信片", flush=True)
    with Pool(10) as p: res = p.map(one, rows, chunksize=4)
    ok  = sum(x for x in res if isinstance(x, int))
    err = [x for x in res if not isinstance(x, int)]
    print(f"新增 ROI {ok} 个;失败 {len(err)}", flush=True)
    for e in err[:10]: print(" ", e)
