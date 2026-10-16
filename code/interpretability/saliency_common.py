# 解释性分析的共用部分(saliency.py、saliency_images.py、saliency_figure.py):训练时预处理与解剖分区。
import os, numpy as np, torch, torch.nn.functional as F
from PIL import Image
KB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ROI = os.environ.get("KNEE_ROI", os.path.join(KB, "data", "knee_roi"))                   # ROI 影像(code/images 流水线的输出)
OUT = os.path.join(KB, "results", "analysis", "saliency")
MODELS = os.environ.get("KSAA_MODELS", os.path.join(KB, "models"))    # 全数据最终模型 m5_full_<backbone>_full.pt
SNAP = os.path.join(KB, "code", "model")        # Mammo-AGE 模型代码
S, R = 448, 112                                                     # 模型输入边长;存档 / 作图用的网格边长
REG = ["Femur, lateral", "Femur, medial", "Joint band, lateral", "Joint band, medial", "Tibia, lateral", "Tibia, medial"]
HALF = 8                                                             # 关节带 = 关节线 ± 8 行(112 网格;约为 ROI 高度的 15%)

def roi_file(seqn, side):                                           # 3 人只有第二张膝片(K2)
    f = os.path.join(ROI, f"{seqn}K1_{side}.png")
    return f if os.path.exists(f) else os.path.join(ROI, f"{seqn}K2_{side}.png")

def pre(f, side):
    """逐项复现 model/knee_v2.py:pct(每图 0.5–99.5 百分位 → [0,1])、左膝水平镜像、双线性缩放到 448。返回 [S,S] 单位灰度张量。"""
    im = np.asarray(Image.open(f)).astype(np.float32)
    lo, hi = np.percentile(im[::2, ::2], [0.5, 99.5]); x = np.clip((im - lo) / max(hi - lo, 1.0), 0.0, 1.0)
    if side == "L": x = x[:, ::-1]
    t = torch.from_numpy(np.ascontiguousarray(x))[None, None]
    if t.shape[-1] != S: t = F.interpolate(t, size=(S, S), mode="bilinear", align_corners=False)
    return t[0, 0]

def pool(t):                                                         # [S,S] → [R,R] 平均池化
    return F.adaptive_avg_pool2d(t[None, None], R)[0, 0]

def joint_line(u):
    """u:[R,R] 单位灰度。关节线 = 中央半宽(列 28–83)行剖面上最深的"谷"(与上下 6 行均值之差最大;行 36–76 内搜索)。
    镜像后外侧(腓骨侧)在图左、内侧在图右。"""
    p = np.convolve(np.asarray(u, np.float32)[:, 28:84].mean(1), np.ones(3) / 3, "same"); r = np.arange(36, 77)
    return int(r[((p[r - 6] + p[r + 6]) / 2 - p[r]).argmax()])

def regions(u):
    """6 个解剖区(股骨 / 关节带 / 胫骨 × 外侧 / 内侧),[6,R,R] 布尔,互不重叠且覆盖全图。"""
    j = joint_line(u); rows = np.arange(R)[:, None]; med = np.arange(R)[None, :] >= R // 2
    band = np.where(rows < j - HALF, 0, np.where(rows <= j + HALF, 1, 2))
    return np.stack([(band == b) & (med == bool(s)) for b in range(3) for s in range(2)]), j

def otsu(u, bins=64):                                               # 骨 / 软组织阈值(单位灰度)
    h, e = np.histogram(np.asarray(u, np.float32), bins=bins, range=(0, 1)); p = h / h.sum(); w = np.cumsum(p)
    mu = np.cumsum(p * (e[:-1] + e[1:]) / 2); s = (mu[-1] * w - mu) ** 2 / (w * (1 - w) + 1e-12)
    return float(e[np.argmax(s) + 1])
