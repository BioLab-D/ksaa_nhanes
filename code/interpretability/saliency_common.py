# Shared code for the interpretability scripts: paths, preprocessing as in training, anatomical regions.
import os, numpy as np, torch, torch.nn.functional as F
from PIL import Image
KB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ROI = os.environ.get("KNEE_ROI", os.path.join(KB, "data", "knee_roi"))                   # ROI images (output of the code/images pipeline)
OUT = os.path.join(KB, "results", "analysis", "saliency")
MODELS = os.environ.get("KSAA_MODELS", os.path.join(KB, "models"))    # final full-data models full_<backbone>.pt
SNAP = os.path.join(KB, "code", "model")        # Mammo-AGE model code
S, R = 448, 112                                                     # model input size; grid size for saved maps and figures
REG = ["Femur, lateral", "Femur, medial", "Joint band, lateral", "Joint band, medial", "Tibia, lateral", "Tibia, medial"]
HALF = 8                                                             # joint band = joint line ± 8 rows (112 grid; about 15% of ROI height)

def roi_file(seqn, side):                                           # 3 participants have only a second knee film (K2)
    f = os.path.join(ROI, f"{seqn}K1_{side}.png")
    return f if os.path.exists(f) else os.path.join(ROI, f"{seqn}K2_{side}.png")

def pre(f, side):
    """Preprocessing of code/model/knee_v2.py (--norm pct): per-image 0.5-99.5 percentiles -> [0,1], left knee mirrored,
    bilinear resize to 448. Returns an [S,S] tensor."""
    im = np.asarray(Image.open(f)).astype(np.float32)
    lo, hi = np.percentile(im[::2, ::2], [0.5, 99.5]); x = np.clip((im - lo) / max(hi - lo, 1.0), 0.0, 1.0)
    if side == "L": x = x[:, ::-1]
    t = torch.from_numpy(np.ascontiguousarray(x))[None, None]
    if t.shape[-1] != S: t = F.interpolate(t, size=(S, S), mode="bilinear", align_corners=False)
    return t[0, 0]

def pool(t):                                                         # [S,S] -> [R,R] average pooling
    return F.adaptive_avg_pool2d(t[None, None], R)[0, 0]

def joint_line(u):
    """Joint line of a [R,R] image: the row (36-76) whose mean intensity over columns 28-83 is furthest below the mean of the
    rows 6 above and 6 below. After mirroring, lateral (fibular side) is on the left and medial on the right."""
    p = np.convolve(np.asarray(u, np.float32)[:, 28:84].mean(1), np.ones(3) / 3, "same")
    r = np.arange(36, 77)
    return int(r[((p[r - 6] + p[r + 6]) / 2 - p[r]).argmax()])

def regions(u):
    """6 non-overlapping regions covering the image (femur / joint band / tibia x lateral / medial), [6,R,R] boolean."""
    j = joint_line(u)
    rows = np.arange(R)[:, None]
    med = np.arange(R)[None, :] >= R // 2
    band = np.where(rows < j - HALF, 0, np.where(rows <= j + HALF, 1, 2))
    return np.stack([(band == b) & (med == bool(s)) for b in range(3) for s in range(2)]), j

def otsu(u, bins=64):                                               # bone / soft-tissue threshold (unit gray)
    h, e = np.histogram(np.asarray(u, np.float32), bins=bins, range=(0, 1))
    p = h / h.sum()
    w = np.cumsum(p)
    mu = np.cumsum(p * (e[:-1] + e[1:]) / 2); s = (mu[-1] * w - mu) ** 2 / (w * (1 - w) + 1e-12)
    return float(e[np.argmax(s) + 1])
