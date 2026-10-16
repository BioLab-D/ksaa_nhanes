"""用启发式自举 YOLO 膝关节检测数据集。

输入 = 整张双膝片(降采样到 1024),输出 = 2 个框(class 0 = 受试者右膝 / 1 = 左膝)。
框由 v4 的关节定位逻辑给出:关节行定位实测 24/24 正确,可直接用作候选;
ROI 边界是启发式的弱项,正是要靠人工审核 + 模型学习来修的部分。

同时输出 suspicion 分,用于把最可能标错的片子排到审核队列前面:
  - 去趋势凹陷强度低 -> 关节定位不可信
  - 左右膝关节行差异大 -> 至少一侧定位错了(站立位双膝关节高度应接近)
  - 腓骨侧别判定与"左半=右膝"的默认假设冲突 -> 可能整片镜像
"""
import os, sys, glob, random, csv, numpy as np
from multiprocessing import Pool
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None

DIR = "/home/data/knee_xray"
OUT = "/home/data/knee_yolo"
SAT = 63000
LONG = 1024
# 框宽(相对整片宽度)的稳健区间,由 300 张实测分布定(中位 0.185);标准化协议片膝宽近乎常数,钳制不损失信号
W_LO, W_HI, W_MED = 0.15, 0.23, 0.185

def _box(x, k):
    k = max(3, int(k) | 1); pad = k // 2
    return np.convolve(np.pad(x, pad, mode="edge"), np.ones(k)/k, mode="valid")

def _longest_run(on, fallback):
    # best 须从零跨度起步,否则更新条件恒为假
    best, cur, s = (0, 0), 0, 0
    for i, v in enumerate(on):
        if v:
            if cur == 0: s = i
            cur += 1
            if cur > best[1]-best[0]: best = (s, i+1)
        else: cur = 0
    return best if best[1] > best[0] else (0, fallback)

def limb_cols(h):
    # 均值定位 + 标准差门控:胶片边缘亮带列均值高但 std 近零,先按 std 排掉,再用均值找最长连续段
    m = np.where(h >= SAT, 0, h).astype(np.float32)
    band = m[int(.25*m.shape[0]):int(.75*m.shape[0]), :]
    pm, ps = band.mean(axis=0), band.std(axis=0)
    on = (pm > pm.max()*0.45) & (ps > ps.max()*0.20)
    return _longest_run(on, h.shape[1])

def knee_cols(h, jr, limb_w):
    b = int(limb_w*0.35)
    m = np.where(h >= SAT, 0, h).astype(np.float32)[max(0,jr-b):min(h.shape[0],jr+b), :]
    pm, ps = m.mean(axis=0), m.std(axis=0)
    r = _longest_run((pm > pm.max()*0.50) & (ps > ps.max()*0.20), h.shape[1])
    return r if r[1]-r[0] >= 50 else (0, h.shape[1])

def joint_row(h, c0, c1):
    m = np.where(h >= SAT, 0, h)[:, c0:c1]
    p = m.mean(axis=1).astype(np.float64); H = len(p)
    dip = _box(p, H*0.10) - _box(p, H*0.008)
    lo, hi = int(.20*H), int(.72*H)
    seg = dip[lo:hi]
    return lo + int(np.argmax(seg)), float(seg.max()/(p.std()+1e-6))

def fibula_side(h, c0, c1, jr):
    w = c1-c0
    r0, r1 = jr+int(w*0.25), min(jr+int(w*0.75), h.shape[0])
    if r1-r0 < 20: return None
    seg = np.where(h >= SAT, 0, h)[r0:r1, c0:c1].mean(axis=0)
    q = max(1, len(seg)//4)
    return "R" if seg[:q].mean() > seg[-q:].mean() else "L"

def window_full(a):
    t = a[(a > 0) & (a < SAT)]
    lo, hi = np.percentile(t if t.size > 1000 else a, [1, 99])
    if hi <= lo: lo, hi = float(a.min()), float(a.max()+1)
    return np.clip((a.astype(np.float32)-lo)/(hi-lo), 0, 1)

def one(p):
    stem = os.path.basename(p)[:-5]
    try:
        a = np.asarray(Image.open(p))
        H, W = a.shape
        prof = np.where(a >= SAT, 0, a).mean(axis=0)
        lo, hi = int(W*.30), int(W*.70)
        mid = lo + int(np.argmin(prof[lo:hi]))

        rows, boxes = [], []
        for cls, (side, off, h) in enumerate((("R", 0, a[:, :mid]), ("L", mid, a[:, mid:]))):
            c0, c1 = limb_cols(h)
            if c1-c0 < 50: c0, c1 = 0, h.shape[1]
            jr, strength = joint_row(h, c0, c1)
            c0, c1 = knee_cols(h, jr, c1-c0)
            fib = fibula_side(h, c0, c1, jr)
            w = c1-c0; cc = (c0+c1)//2 + off
            boxes.append([cls, cc, jr, w*1.16])       # 尺寸稍后统一钳制
            rows.append(dict(side=side, jr=jr, w=w, strength=strength, fib=fib))

        # 中心可靠、尺寸不可靠:左右膝取同一宽度并钳到稳健区间
        bw = float(np.clip(np.mean([b[3] for b in boxes])/W, W_LO, W_HI)) * W
        for b in boxes: b[3] = bw

        dj = abs(rows[0]["jr"] - rows[1]["jr"])
        lat_ok = (rows[0]["fib"] == "R" and rows[1]["fib"] == "L")
        # 侧别不进评分:左半 = 受试者右膝(审核图已证实);腓骨判据无信息量
        susp = (0.0 if rows[0]["strength"] > .35 else 1.0) + \
               (0.0 if rows[1]["strength"] > .35 else 1.0) + \
               (1.5 if dj > 250 else 0.0)

        sc = LONG / max(H, W)
        img = Image.fromarray((window_full(a)*255).astype(np.uint8)) \
                   .resize((int(W*sc), int(H*sc)), Image.LANCZOS)
        img.save(f"{OUT}/images/{stem}.png")
        with open(f"{OUT}/labels/{stem}.txt", "w") as f:
            for cls, cx, cy, side_px in boxes:
                f.write(f"{cls} {cx/W:.6f} {cy/H:.6f} {side_px/W:.6f} {side_px/H:.6f}\n")
        return dict(stem=stem, ok=1, mid=mid, W=W, H=H, dj=dj,
                    s_r=round(rows[0]["strength"],3), s_l=round(rows[1]["strength"],3),
                    fib_r=rows[0]["fib"], fib_l=rows[1]["fib"],
                    lat_ok=int(lat_ok), susp=round(susp,2))
    except Exception as e:
        return dict(stem=stem, ok=0, err=f"{type(e).__name__}: {e}", susp=99)

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    for d in ("images", "labels", "review"):
        os.makedirs(f"{OUT}/{d}", exist_ok=True)
    files = sorted(glob.glob(DIR+"/*.tiff"))
    random.seed(7)
    samp = random.sample(files, min(n, len(files)))
    with Pool(14) as pool:
        res = pool.map(one, samp, chunksize=4)
    keys = ["stem","ok","W","H","mid","dj","s_r","s_l","fib_r","fib_l","lat_ok","susp","err"]
    with open(f"{OUT}/meta.csv","w",newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); wtr.writeheader()
        for r in res: wtr.writerow(r)
    ok = [r for r in res if r["ok"]]
    print(f"成功 {len(ok)} / {len(res)}")
    if ok:
        import statistics as st
        print(f"  可疑分>0 的片子: {sum(1 for r in ok if r['susp']>0)}")
        print(f"  侧别判定与默认冲突: {sum(1 for r in ok if not r['lat_ok'])}")
        print(f"  左右关节行差 中位 {st.median(r['dj'] for r in ok):.0f}  "
              f"90分位 {sorted(r['dj'] for r in ok)[int(.9*len(ok))]:.0f}")
