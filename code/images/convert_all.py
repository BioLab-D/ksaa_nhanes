"""把全部 4,908 张膝片转成 1024px 8-bit PNG,作为检测器输入。

窗宽只用组织像素统计:直曝背景和铅标记已饱和到 65535,
纳入分位数会把骨与软组织压到很窄的灰阶里。
已存在的跳过(前面 800 张候选已生成)。
"""
import os, glob, sys, numpy as np
from multiprocessing import Pool
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
DIR, OUT, SAT, LONG = "/home/data/knee_xray", "/home/data/knee_yolo/images", 63000, 1024

def one(p):
    stem = os.path.basename(p)[:-5]
    dst = f"{OUT}/{stem}.png"
    if os.path.exists(dst): return "skip"
    try:
        a = np.asarray(Image.open(p)); H, W = a.shape
        t = a[(a > 0) & (a < SAT)]
        lo, hi = np.percentile(t if t.size > 1000 else a, [1, 99])
        if hi <= lo: lo, hi = float(a.min()), float(a.max()+1)
        b = np.clip((a.astype(np.float32)-lo)/(hi-lo), 0, 1)
        sc = LONG/max(H, W)
        Image.fromarray((b*255).astype(np.uint8)) \
             .resize((int(W*sc), int(H*sc)), Image.LANCZOS).save(dst)
        return "ok"
    except Exception as e:
        return f"err {stem} {type(e).__name__}"

if __name__ == "__main__":
    files = sorted(glob.glob(DIR+"/*.tiff"))
    print(f"总片数 {len(files)}  已存在 {len(glob.glob(OUT+'/*.png'))}", flush=True)
    with Pool(12) as pool:
        r = pool.map(one, files, chunksize=8)
    from collections import Counter
    c = Counter(x if x in ("ok","skip") else "err" for x in r)
    print(dict(c), flush=True)
    for x in r:
        if x not in ("ok","skip"): print(" ", x, flush=True)
    print(f"最终 PNG 数: {len(glob.glob(OUT+'/*.png'))}", flush=True)
