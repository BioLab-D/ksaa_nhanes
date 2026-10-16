"""按检测框从原始 16-bit TIFF 裁高分辨率膝关节 ROI。

存 16-bit 而非窗宽后的 8-bit:重裁一遍要把 4,908 张 38MB 的 TIFF 全部重读
(约 1 小时),而磁盘有 2.9TB 余量。保留原始动态范围,把"用什么窗宽"这个
下游还会反复调的决定留到训练时再做,不锁死在预处理里。

框外扩 1.15 倍留少量解剖上下文;超出片边的用 edge 填充而非补零,
避免在 ROI 边缘制造人工的强梯度。
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
    print(f"待裁 {len(stems)} 张片 / {sum(len(v) for v in boxes.values())} 个膝", flush=True)
    with Pool(12) as pool:
        res = pool.map(one, stems, chunksize=8)
    flat = [x for sub in res for x in sub]
    from collections import Counter
    print(dict(Counter(x if x in ("ok","skip") else "err" for x in flat)))
    for x in flat:
        if x not in ("ok","skip"): print(" ", x)
